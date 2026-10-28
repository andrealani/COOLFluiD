#!/usr/bin/env bash
# Initialize large COCONUT regression test datasets.
# Run from the COOLFluiD repository root.

# Reusable for any testcase and Zenodo record.
fetch_and_extract_from_zenodo() {
    local dataset="$1" url="$2"
    local archive="${dataset}.tar.xz" part="${dataset}.tar.xz.part"
    local marker=".${dataset}.extracted" staging=".${dataset}.extracting"
    local attempt max_attempts=5 valid=0
    local -a wget_progress=(--quiet)
    if [[ -t 2 ]]; then
        wget_progress=(--progress=bar:force:noscroll)
    fi

    if [[ -f "$marker" && -d "$dataset" ]]; then
        echo "$dataset already extracted; skipping."
        return 0
    fi
    echo "Preparing $dataset..."

    # A directory without its marker may be the result of an interrupted extraction.
    # Don't overwrite unknown data automatically.
    if [[ -d "$dataset" ]]; then
        echo "ERROR: $dataset exists without $marker; inspect it before retrying." >&2
        return 1
    fi
    # This staging directory belongs exclusively to this function.
    rm -rf -- "$staging"

    for ((attempt=1; attempt<=max_attempts; attempt++)); do
        if [[ -f "$archive" ]]; then
            if xz -t "$archive" && tar -tJf "$archive" >/dev/null; then
                valid=1
                break
            fi
            echo "Invalid $archive; deleting it and downloading again." >&2
            rm -f -- "$archive"
        fi

        echo "Downloading $dataset (attempt $attempt/$max_attempts)..."
        # --continue resumes a partial download; --tries handles transient failures.
        # Show a compact bar on terminals; stay quiet in logs.
        if wget "${wget_progress[@]}" --continue --tries=3 --timeout=30 \
             -O "$part" "$url"; then
            if xz -t "$part" && tar -tJf "$part" >/dev/null; then
                mv -f -- "$part" "$archive"
                valid=1
                break
            fi
            echo "Downloaded archive failed validation; restarting from zero." >&2
            rm -f -- "$part"
        else
            echo "Download interrupted; retaining $part for the next attempt." >&2
        fi
    done
    if (( ! valid )); then
        echo "ERROR: Could not obtain a valid $archive after $max_attempts attempts." >&2
        return 1
    fi

    if ! mkdir -p -- "$staging"; then
        echo "ERROR: Could not create staging directory $staging." >&2
        return 1
    fi
    # Approximate extraction progress, using standard tools only (no pv).
    # xz --robot --list: totals STREAMS BLOCKS COMPRESSED UNCOMPRESSED ...
    local total_bytes total_kb monitor_pid="" extraction_ok=0
    total_bytes=$(xz --robot --list "$archive" 2>/dev/null | awk '$1 == "totals" {print $5}')
    # In xz robot output, fields are: totals, streams, blocks, compressed, uncompressed.
    # The fifth field is the uncompressed byte count.
    if [[ "$total_bytes" =~ ^[0-9]+$ ]] && (( total_bytes > 0 )); then
        total_kb=$(( (total_bytes + 1023) / 1024 ))
    else
        total_kb=0
    fi
    (
        used_kb=0
        tick=0
        while :; do
            used_kb=$(du -sk "$staging" 2>/dev/null | awk '{print $1}')
            [[ "$used_kb" =~ ^[0-9]+$ ]] || used_kb=0
            if (( total_kb > 0 )); then
                pct=$(( used_kb * 100 / total_kb ))
                (( pct > 99 )) && pct=99
                filled=$(( pct / 5 ))
                printf -v bar '%*s' "$filled" ''
                bar=${bar// /#}
                printf '\rExtracting %s: [%-20s] %2d%% (approx.)' "$dataset" "$bar" "$pct" >&2
            else
                printf '\rExtracting %s: %s KB written (elapsed %ds)' "$dataset" "$used_kb" "$tick" >&2
            fi
            tick=$((tick + 1))
            sleep 1
        done
    ) &
    monitor_pid=$!

    if tar -xJf "$archive" -C "$staging"; then
        extraction_ok=1
    fi
    kill "$monitor_pid" 2>/dev/null || :
    wait "$monitor_pid" 2>/dev/null || :
    if (( extraction_ok )); then
        printf '\rExtracting %s: [####################] 100%% (complete)\n' "$dataset" >&2
    else
        printf '\n' >&2
    fi
    if (( ! extraction_ok )); then
        echo "ERROR: extraction of $archive failed." >&2
        rm -rf -- "$staging"
        return 1
    fi
    if [[ ! -d "$staging/$dataset" ]]; then
        echo "ERROR: archive does not contain the expected $dataset/ directory." >&2
        rm -rf -- "$staging"
        return 1
    fi
    if ! mv -- "$staging/$dataset" "$dataset"; then
        echo "ERROR: Could not move $dataset into place." >&2
        return 1
    fi
    if ! touch -- "$marker"; then
        echo "ERROR: Could not create completion marker $marker." >&2
        return 1
    fi
    rm -rf -- "$staging"
    if ! rm -f -- "$archive"; then
        echo "WARNING: Could not delete $archive." >&2
    fi
    echo "$dataset extraction complete."
}

# Unpack/fetch large data (e.g. meshes, maps) for the corresponding testcases 

set -e

if [ -f plugins/MHD/testcases/COCONUT/Dipole_WTD/30Rs.CFmesh.xz ]; then
    unxz plugins/MHD/testcases/COCONUT/Dipole_WTD/30Rs.CFmesh.xz
fi
if [ -f plugins/MHD/testcases/COCONUT/Dipole_Shifted/2008_lmax25.dat.xz ]; then
    unxz plugins/MHD/testcases/COCONUT/Dipole_Shifted/2008_lmax25.dat.xz
fi
if [ -f plugins/MHD/testcases/COCONUT/Dipole_Shifted/2Rslvl6.CFmesh.xz ]; then
    unxz plugins/MHD/testcases/COCONUT/Dipole_Shifted/2Rslvl6.CFmesh.xz
fi

(
    cd plugins/MHD/testcases/COCONUT/Eclipse/Mesh
    if compgen -G 'corona_720.CFmesh.*.xz' > /dev/null; then
        unxz corona_720.CFmesh.*.xz
        cat corona_720.CFmesh.?? > corona_720.CFmesh
    fi
)
(
    cd plugins/MHD/testcases/COCONUT/Eclipse/MapData
    if compgen -G '*.xz' > /dev/null; then
        unxz map_gong_lmax25*.xz
    fi
)
(
    cd plugins/MHD/testcases/COCONUT/Eclipse600/Mesh
    if compgen -G 'L6corona_PPDecEBC2_2024zqsLax025NHom.CFmesh.*.xz' > /dev/null; then
        unxz L6corona_PPDecEBC2_2024zqsLax025NHom.CFmesh.*.xz
        cat L6corona_PPDecEBC2_2024zqsLax025NHom.CFmesh.?? > L6corona_PPDecEBC2_2024zqsLax025NHom.CFmesh
    fi
)
(
    cd plugins/MHD/testcases/COCONUT/Eclipse600/MapData/First600
    if compgen -G '*.xz' > /dev/null; then
        unxz map_gong_lmax50_20240*.xz
    fi
)
(
    cd plugins/MHD/testcases/COCONUT/SunToEarthCR2282
    zenodo_record='https://zenodo.org/records/23600389/files'
    printf '\rplugins/MHD/testcases/COCONUT/SunToEarthCR2282 downloading from:\t%s\n' "[ $zenodo_record ]" >&2
    fetch_and_extract_from_zenodo Mesh "${zenodo_record}/Mesh.tar.xz?download=1"
    fetch_and_extract_from_zenodo MapData "${zenodo_record}/MapData.tar.xz?download=1"
)

if [ -f plugins/RadiativeTransfer/testcases/SolarCorona/corona_fullMHD.CFmesh.xz ]; then
    unxz plugins/RadiativeTransfer/testcases/SolarCorona/corona_fullMHD.CFmesh.xz
fi

echo '######## init_testcases done! ########'
