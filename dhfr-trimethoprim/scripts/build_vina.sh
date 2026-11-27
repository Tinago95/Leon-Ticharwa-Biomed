#!/bin/bash
# Build AutoDock Vina 1.2.7 from source on macOS (Apple Silicon) and install to ~/.local/bin.
#
# Why from source: there is no Homebrew formula, and `pip install vina` fails because its
# setup.py searches only "conda env", /usr/local/include and /usr/include - it ignores
# BOOST_ROOT/CPATH/CPPFLAGS, so it never finds Homebrew Boost at /opt/homebrew.
#
# Do not substitute the prebuilt 1.1.2 macOS binary: it is unsigned, x86_64, and macOS
# deletes it on execution. It is also the wrong version (the project plan targets 1.2).
#
# Usage:  bash scripts/build_vina.sh

set -euo pipefail

VERSION="1.2.7"
PREFIX="${HOME}/.local/bin"
WORKDIR="$(mktemp -d)"
trap 'rm -rf "${WORKDIR}"' EXIT

command -v brew >/dev/null || { echo "Homebrew required"; exit 1; }

if ! brew list --versions boost >/dev/null 2>&1; then
    echo "==> installing boost"
    brew install boost
fi
echo "==> boost $(brew list --versions boost)"

echo "==> fetching Vina ${VERSION}"
cd "${WORKDIR}"
curl -fsSL -o "vina-${VERSION}.tar.gz" \
    "https://github.com/ccsb-scripps/AutoDock-Vina/archive/refs/tags/v${VERSION}.tar.gz"
tar xzf "vina-${VERSION}.tar.gz"

echo "==> building (Makefile resolves brew --prefix itself)"
cd "AutoDock-Vina-${VERSION}/build/mac/release"
make -j"$(sysctl -n hw.ncpu)"

echo "==> installing to ${PREFIX}"
mkdir -p "${PREFIX}"
install -m 755 vina vina_split "${PREFIX}/"

echo "==> verifying"
file "${PREFIX}/vina"
"${PREFIX}/vina" --version

case ":${PATH}:" in
    *":${PREFIX}:"*) ;;
    *) echo "NOTE: ${PREFIX} is not on PATH - add it to ~/.zshrc" ;;
esac
