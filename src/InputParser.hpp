// ===========================================================================
//  InputParser.hpp — Axom Inlet-based YAML input parser
// ===========================================================================
#pragma once

#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "Functions.hpp"
#include "InputData.hpp"

// ===========================================================================
//  FromInlet specializations  (must live in the global namespace)
//
//  NOTE on Axom Inlet collection retrieval
//  
//  Container::get<T>() for collection types (array/dict/vector) takes *no*
//  name argument.  To read a named sub-collection use the Proxy path:
//
//    base["field"].get<std::vector<double>>()
//
//  Container::isUserProvided(name) correctly returns false when the user
//  omitted an optional array, even if the field was declared in the schema.
// ===========================================================================

// ---------------------------------------------------------------------------
//  ShellInput
// ---------------------------------------------------------------------------
template <>
struct FromInlet<ShellInput>
{
    ShellInput operator()(const axom::inlet::Container& base) const
    {
        ShellInput s;
        s.id           = base["id"];
        s.outer_radius = base["outer_radius"];
        s.num_cells    = base["num_cells"];
        s.material     = static_cast<std::string>(base["material"]);
        s.source       = static_cast<std::string>(base["source"]);
        return s;
    }
};

// ---------------------------------------------------------------------------
//  MaterialInput
// ---------------------------------------------------------------------------
template <>
struct FromInlet<MaterialInput>
{
    MaterialInput operator()(const axom::inlet::Container& base) const
    {
        MaterialInput m;
        m.id              = static_cast<std::string>(base["id"]);
        m.sigma_t         = base["sigma_t"];
        m.sigma_a         = base["sigma_a"];
        m.scattering_order = base["scattering_order"];  // has defaultValue(0)

        if (m.scattering_order > 0 && base.isUserProvided("sigma_s_moments"))
        {
            m.sigma_s_moments =
                base["sigma_s_moments"].get<std::vector<double>>();
        }
        return m;
    }
};

// ---------------------------------------------------------------------------
//  DistributedSourceInput
// ---------------------------------------------------------------------------
template <>
struct FromInlet<DistributedSourceInput>
{
    DistributedSourceInput operator()(const axom::inlet::Container& base) const
    {
        DistributedSourceInput d;
        d.id        = static_cast<std::string>(base["id"]);
        d.normalize = base["normalize"];   // has defaultValue(false)

        if (base.isUserProvided("strength"))
            d.strength = base["strength"];

        if (base.isUserProvided("per_cell_strength"))
            d.per_cell_strength =
                base["per_cell_strength"].get<std::vector<double>>();

        return d;
    }
};

// ===========================================================================
//  InputParser
// ===========================================================================
class InputParser
{
public:
    explicit InputParser(const std::string& filepath) : m_filepath(filepath) {}

    // -----------------------------------------------------------------------
    //  parse() — open file, define schema, verify, extract → InputData
    // -----------------------------------------------------------------------
    InputData parse()
    {
        namespace inlet = axom::inlet;

        auto reader = std::make_unique<inlet::YAMLReader>();
        if (!reader->parseFile(m_filepath))
            throw std::runtime_error(
                "InputParser: failed to read YAML file: " + m_filepath);

        inlet::Inlet inlet_obj(std::move(reader), /*docEnabled=*/false);
        define_schema(inlet_obj);

        std::vector<inlet::VerificationError> errors;
        if (!inlet_obj.verify(&errors))
        {
            std::ostringstream oss;
            oss << "InputParser: verification failed for '"
                << m_filepath << "':\n";
            for (const auto& e : errors)
                oss << "  " << e.message << "\n";
            throw std::runtime_error(oss.str());
        }

        return extract(inlet_obj);
    }

private:
    // -----------------------------------------------------------------------
    //  define_schema - register every field Inlet should look for
    // -----------------------------------------------------------------------
    static void define_schema(axom::inlet::Inlet& inlet)
    {
        //  name
        inlet.addString("name", "Name of the problem").defaultValue("");

        //  geometry 
        auto& geom   = inlet.addStruct("geometry", "Problem geometry");
        auto& shells = geom.addStructArray("shells",
                           "Radial shells listed innermost -> outermost");
        shells.addInt   ("id",           "Shell identifier").isRequired();
        shells.addDouble("outer_radius", "Outer radius [cm]").isRequired();
        shells.addInt   ("num_cells",    "Number of spatial cells").isRequired();
        shells.addString("material",     "Material id").isRequired();
        shells.addString("source",       "Distributed source id, or 'none'");

        //  quadrature 
        auto& quad = inlet.addStruct("angular_quadrature", "Angular quadrature");
        quad.addInt("sn_order", "Sn order N (must be even)").isRequired();

        //  materials 
        auto& mats = inlet.addStructArray("materials",
                         "Material cross-section data");
        mats.addString("id",              "Material identifier").isRequired();
        mats.addDouble("sigma_t",         "Total cross section [1/cm]").isRequired();
        mats.addDouble("sigma_a",         "Absorption cross section [1/cm]").isRequired();
        mats.addInt   ("scattering_order","Pn scattering order K (0=isotropic)")
                      .defaultValue(0);
        mats.addDoubleArray("sigma_s_moments",
                            "Scattering moments sigma_sk for k=1..K");

        //  sources 
        auto& sources = inlet.addStruct("sources", "Neutron sources");

        auto& dist = sources.addStructArray("distributed",
                               "Distributed volumetric sources");
        dist.addString("id",    "Source identifier").isRequired();
        dist.addDouble("strength",
                       "Uniform source strength [p/cm3/s]");
        dist.addDoubleArray("per_cell_strength",
                            "Per-cell strengths (overrides 'strength')");
        dist.addBool("normalize",
                     "Normalize to unit volume integral").defaultValue(false);

        auto& bnd = sources.addStruct("boundary",
                        "Outer-boundary incident angular flux");
        bnd.addDouble("isotropic_flux",
                      "Isotropic incident angular flux ψ (applied to all incoming directions)")
                     .defaultValue(0.0);
        bnd.addDouble("scalar_flux",
                      "Isotropic incident scalar flux φ; solver converts to ψ = φ/2")
                     .defaultValue(0.0);
        bnd.addDoubleArray("per_direction",
                           "Per-direction psi_m for mu_m < 0 (length = N/2)");
        bnd.addBool("normalize",
                    "Normalize to unit half-range current").defaultValue(false);

        //  solver 
        auto& solver = inlet.addStruct("solver",
                           "Source-iteration solver settings");
        solver.addString("acceleration",
                         "Acceleration: none | dsa").defaultValue("none");
        solver.addDouble("convergence_tolerance",
                         "Scalar-flux L-inf convergence tolerance")
                        .defaultValue(1.0e-4);
        solver.addInt("max_iterations",
                      "Maximum source iterations").defaultValue(500);

        //  output 
        auto& output = inlet.addStruct("output", "Output controls");
        output.addString("directory",
                         "Results directory").defaultValue("./results/");
        output.addBool("scalar_flux_csv",
                       "Write cell-centered scalar flux into csv").defaultValue(true);
        output.addBool("scalar_flux_pdv",
                       "Write cell-centered scalar flux into file for pdv").defaultValue(true);
        output.addBool("balance_table",
                       "Write global balance table").defaultValue(true);
        output.addBool("angular_flux_boundary",
                       "Write angular flux at outer boundary").defaultValue(true);
        output.addBool("starting_direction_origin",
                       "Write psi_{1/2} at r=0").defaultValue(true);
    }

    // -----------------------------------------------------------------------
    //  extract — pull parsed values out of Inlet into InputData
    // -----------------------------------------------------------------------
    static InputData extract(axom::inlet::Inlet& inlet)
    {
        InputData data;
        
        // name
        data.name = inlet.get<std::string>("name");

        //  geometry/shells
        data.shells =
            inlet["geometry/shells"].get<std::vector<ShellInput>>();
        std::sort(data.shells.begin(), data.shells.end(),
                  [](const ShellInput& a, const ShellInput& b)
                  { return a.id < b.id; });

        //  quadrature 
        data.quadrature.order = inlet.get<int>("angular_quadrature/sn_order");
        GaussLegendreQuadrature(data.quadrature.order,
                                data.quadrature.mu,
                                data.quadrature.w);

        //  materials 
        data.materials =
            inlet["materials"].get<std::vector<MaterialInput>>();

        //  sources/distributed (entire section is optional) 
        if (inlet.isUserProvided("sources/distributed"))
        {
            data.distributed_sources =
                inlet["sources/distributed"]
                    .get<std::vector<DistributedSourceInput>>();
        }

        //  sources/boundary
        // Scalar fields have defaultValue() so always safe to read.
        data.boundary.isotropic_flux =
            inlet.get<double>("sources/boundary/isotropic_flux");
        data.boundary.scalar_flux =
            inlet.get<double>("sources/boundary/scalar_flux");
        data.boundary.normalize =
            inlet.get<bool>("sources/boundary/normalize");

        if (inlet.isUserProvided("sources/boundary/per_direction"))
        {
            data.boundary.per_direction =
                inlet["sources/boundary/per_direction"]
                    .get<std::vector<double>>();
        }

        //  solver 
        data.solver.acceleration =
            inlet.get<std::string>("solver/acceleration");
        data.solver.convergence_tolerance =
            inlet.get<double>("solver/convergence_tolerance");
        data.solver.max_iterations =
            inlet.get<int>("solver/max_iterations");

        //  output 
        data.output.directory =
            inlet.get<std::string>("output/directory");
        data.output.scalar_flux_csv =
            inlet.get<bool>("output/scalar_flux_csv");
        data.output.scalar_flux_pdv =
            inlet.get<bool>("output/scalar_flux_pdv");
        data.output.balance_table =
            inlet.get<bool>("output/balance_table");
        data.output.angular_flux_boundary =
            inlet.get<bool>("output/angular_flux_boundary");
        data.output.starting_direction_origin =
            inlet.get<bool>("output/starting_direction_origin");

        return data;
    }

    std::string m_filepath;
};
