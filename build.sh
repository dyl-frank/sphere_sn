#!/usr/bin/env bash
# ===========================================================================
#  build.sh — Automated build script for 1D Spherical Sn Transport Code
#
#  Usage:
#    ./build.sh              # full build (install deps if needed + compile)
#    ./build.sh --deps-only  # only install/load dependencies
#    ./build.sh --no-deps    # skip dependency checks, just cmake + make
#    ./build.sh --clean      # wipe build directory and rebuild from scratch
#    ./build.sh --help       # print usage
# ===========================================================================

set -euo pipefail

# ── Color output helpers ──────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

info()    { echo -e "${CYAN}[INFO]${NC}    $*"; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $*"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}    $*"; }
error()   { echo -e "${RED}[ERROR]${NC}   $*"; exit 1; }
header()  { echo -e "\n${BOLD}═══════════════════════════════════════════════════${NC}";
            echo -e "${BOLD}  $*${NC}";
            echo -e "${BOLD}═══════════════════════════════════════════════════${NC}\n"; }

# ── Project paths ─────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${SCRIPT_DIR}"
BUILD_DIR="${PROJECT_DIR}/build"
SPACK_ROOT="${SPACK_ROOT:-${HOME}/spack}"

# ── Default settings ──────────────────────────────────────────────────────
BUILD_TYPE="${BUILD_TYPE:-Release}"
NUM_JOBS="${NUM_JOBS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}"
INSTALL_DEPS=true
DO_BUILD=true
CLEAN_BUILD=false

# ── Parse command-line flags ──────────────────────────────────────────────
print_usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Options:
  --deps-only     Only install and load dependencies (skip build)
  --no-deps       Skip dependency checks (assume already loaded)
  --clean         Remove build directory and rebuild from scratch
  --debug         Build in Debug mode instead of Release
  --jobs N        Number of parallel make jobs (default: ${NUM_JOBS})
  --help          Show this help message

Environment variables:
  SPACK_ROOT      Path to Spack installation (default: ~/spack)
  BUILD_TYPE      CMake build type (default: Release)
  NUM_JOBS        Parallel jobs for make (default: nproc)
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --deps-only)  INSTALL_DEPS=true;  DO_BUILD=false;  shift ;;
        --no-deps)    INSTALL_DEPS=false; DO_BUILD=true;   shift ;;
        --clean)      CLEAN_BUILD=true;   shift ;;
        --debug)      BUILD_TYPE="Debug"; shift ;;
        --jobs)       NUM_JOBS="$2";      shift 2 ;;
        --help|-h)    print_usage ;;
        *)            error "Unknown option: $1\n  Run with --help for usage." ;;
    esac
done

# ══════════════════════════════════════════════════════════════════════════
#  Step 1: Ensure Spack is available
# ══════════════════════════════════════════════════════════════════════════
setup_spack() {
    header "Step 1: Spack"

    # Check if spack command is already available
    if command -v spack &>/dev/null; then
        info "Spack already available: $(spack --version)"
        return 0
    fi

    # Try sourcing from SPACK_ROOT
    if [[ -f "${SPACK_ROOT}/share/spack/setup-env.sh" ]]; then
        info "Sourcing Spack from ${SPACK_ROOT}"
        source "${SPACK_ROOT}/share/spack/setup-env.sh"
        return 0
    fi

    # Spack not found — offer to install
    warn "Spack not found at ${SPACK_ROOT}"
    read -rp "Install Spack to ${SPACK_ROOT}? [Y/n] " response
    response="${response:-Y}"
    if [[ ! "${response}" =~ ^[Yy] ]]; then
        error "Spack is required. Install it or set SPACK_ROOT and retry."
    fi

    info "Cloning Spack..."
    git clone -c feature.manyFiles=true \
        https://github.com/spack/spack.git "${SPACK_ROOT}"

    source "${SPACK_ROOT}/share/spack/setup-env.sh"

    info "Detecting compilers..."
    spack compiler find

    success "Spack installed and configured."
    echo ""
    warn "To make Spack available in future shells, add this to your ~/.bashrc:"
    echo "    source ${SPACK_ROOT}/share/spack/setup-env.sh"
    echo ""
}

# ══════════════════════════════════════════════════════════════════════════
#  Step 2: Install dependencies via Spack
# ══════════════════════════════════════════════════════════════════════════
install_deps() {
    header "Step 2: Dependencies"

    # ── Eigen ─────────────────────────────────────────────────────────────
    if spack find --no-groups eigen 2>/dev/null | grep -q eigen; then
        info "Eigen already installed."
    else
        info "Installing Eigen (header-only — this is fast)..."
        spack install eigen
        success "Eigen installed."
    fi

    # ── Axom ──────────────────────────────────────────────────────────────
    if spack find --no-groups axom 2>/dev/null | grep -q axom; then
        info "Axom already installed."
    else
        info "Installing Axom (this may take 15-30 minutes)..."
        info "Building with: ~mpi ~fortran ~examples ~tests (trimmed for faster build)"
        spack install axom~mpi~fortran~examples~tests
        success "Axom installed."
    fi

    # ── Load both into the current environment ────────────────────────────
    info "Loading packages into environment..."

    # EIGEN3_DIR="$(spack location -i eigen)/share/eigen3/cmake"
    # AXOM_DIR="$(spack location -i axom)/lib/cmake/axom"
    success "All dependencies installed and loaded."
}

# ══════════════════════════════════════════════════════════════════════════
#  Step 3: Locate dependency paths for CMake
# ══════════════════════════════════════════════════════════════════════════
find_dependency_paths() {
    header "Step 3: Resolving paths"

    # ── Eigen ─────────────────────────────────────────────────────────────
    EIGEN_PREFIX="$(spack location -i eigen)"
    EIGEN3_DIR="${EIGEN_PREFIX}/share/eigen3/cmake"

    # Fallback: search for the config file
    if [[ ! -f "${EIGEN3_DIR}/Eigen3Config.cmake" ]]; then
        info "Searching for Eigen3Config.cmake..."
        EIGEN3_DIR="$(find "${EIGEN_PREFIX}" -name 'Eigen3Config.cmake' -printf '%h' -quit 2>/dev/null || true)"
        [[ -z "${EIGEN3_DIR}" ]] && error "Could not find Eigen3Config.cmake under ${EIGEN_PREFIX}"
    fi
    info "Eigen3_DIR:  ${EIGEN3_DIR}"

    # ── Axom ──────────────────────────────────────────────────────────────
    AXOM_PREFIX="$(spack location -i axom)"
    AXOM_DIR="${AXOM_PREFIX}/lib/cmake/axom"

    # Fallback: search for the config file
    if [[ ! -f "${AXOM_DIR}/axom-config.cmake" ]] && [[ ! -f "${AXOM_DIR}/axomConfig.cmake" ]]; then
        info "Searching for axom config cmake file..."
        AXOM_DIR="$(find "${AXOM_PREFIX}" -name 'axom*onfig.cmake' -printf '%h' -quit 2>/dev/null || true)"
        [[ -z "${AXOM_DIR}" ]] && error "Could not find axomConfig.cmake under ${AXOM_PREFIX}"
    fi
    info "axom_DIR:    ${AXOM_DIR}"

    success "Dependency paths resolved."
    
    # ── MFEM ──────────────────────────────────────────────────────────────
    MFEM_PREFIX="$(spack location -i mfem)"
    MFEM_CONFIG="${MFEM_PREFIX}/share/mfem/config.mk"

    # Verify config.mk exists
    if [[ ! -f "${MFEM_CONFIG}" ]]; then
        info "Searching for MFEM config.mk..."
        MFEM_CONFIG="$(find "${MFEM_PREFIX}" -name 'config.mk' -path '*mfem*' -print -quit 2>/dev/null || true)"
        [[ -z "${MFEM_CONFIG}" ]] && error "Could not find MFEM config.mk under ${MFEM_PREFIX}"
    fi

    info "MFEM prefix: ${MFEM_PREFIX}"
    info "MFEM config: ${MFEM_CONFIG}"
}

# ══════════════════════════════════════════════════════════════════════════
#  Step 4: Configure and build with CMake
# ══════════════════════════════════════════════════════════════════════════
build_project() {
    header "Step 4: Build (${BUILD_TYPE}, ${NUM_JOBS} jobs)"

    # ── Clean if requested ────────────────────────────────────────────────
    if [[ "${CLEAN_BUILD}" == true ]] && [[ -d "${BUILD_DIR}" ]]; then
        info "Removing existing build directory..."
        rm -rf "${BUILD_DIR}"
    fi

    mkdir -p "${BUILD_DIR}"
    cd "${BUILD_DIR}"

    # ── CMake configure ───────────────────────────────────────────────────
    info "Running CMake configure..."
    cmake "${PROJECT_DIR}" \
        -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
        -DEigen3_DIR="${EIGEN3_DIR}" \
        -Daxom_DIR="${AXOM_DIR}" \
        -DMFEM_DIR="${MFEM_PREFIX}"
    success "CMake configuration complete."

    # ── Build ─────────────────────────────────────────────────────────────
    info "Compiling with ${NUM_JOBS} parallel jobs..."
    cmake --build . --parallel "${NUM_JOBS}"

    success "Build complete."

    # ── Summary ───────────────────────────────────────────────────────────
    echo ""
    info "Executable: ${BUILD_DIR}/sphere_sn"
    echo ""
}

# ══════════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════════
main() {
    header "1D Spherical Sn Transport — Build Script"

    info "Project:    ${PROJECT_DIR}"
    info "Build dir:  ${BUILD_DIR}"
    info "Build type: ${BUILD_TYPE}"
    info "Jobs:       ${NUM_JOBS}"
    echo ""

    if [[ "${INSTALL_DEPS}" == true ]]; then
        setup_spack
        install_deps
    fi

    if [[ "${DO_BUILD}" == true ]]; then
        find_dependency_paths
        build_project
    fi

    header "Done"
}

main
