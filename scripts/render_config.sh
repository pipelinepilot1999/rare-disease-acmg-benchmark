#!/usr/bin/env bash
# Render the Exomiser property files for this machine.
#
# Exomiser is Java and will not expand environment variables inside a properties
# file, so the data directory has to be written in literally. The templates are
# tracked; the rendered copies are not, because they are machine-specific.
set -euo pipefail
HOME_DIR="${RDACMG_HOME:-$HOME/rdacmg}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
for t in "$REPO"/config/*.properties.template; do
  out="${t%.template}"
  sed "s|__RDACMG_HOME__|$HOME_DIR|g" "$t" | grep -v '^# Rendered by\|^# from \$RDACMG_HOME' > "$out"
  echo "rendered $(basename "$out")  (data-directory: $HOME_DIR/resources/exomiser/data)"
done
