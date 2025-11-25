#!/usr/bin/env bash
# Usage: ./time_logs_to_tsv_fast.sh /path/to/timelogs > combined.tsv

logdir="$1"

# Print header
echo -e "file\ttool\tgene\tuser_time\tsystem_time\tcpu_percent\twall_time\twall_time_sec\tmax_rss_MB"

# Function: parse one file with awk (faster than greps)
parse_file() {
    f="$1"
    base=$(basename "$f" .time)
    tool="${base%%_*}"
    gene="${base#*_}"

    awk -v file="$base.time" -v tool="$tool" -v gene="$gene" '
        function to_seconds(t,   n,a) {
            n = split(t, a, ":")
            if (n == 3) return a[1]*3600 + a[2]*60 + a[3]
            else if (n == 2) return a[1]*60 + a[2]
            else return t
        }
        /User time \(seconds\):/    {user=$4}
        /System time \(seconds\):/  {sys=$4}
        /Percent of CPU this job got:/ {cpu=$7; sub("%","",cpu)}
        /Elapsed \(wall clock\) time/ {wall=$NF; wall_sec=to_seconds(wall)}
        /Maximum resident set size \(kbytes\):/ {mem=$6; mem_MB=sprintf("%.2f",mem/1024)}
        END {
            if (user=="") user="NA"
            if (sys=="") sys="NA"
            if (cpu=="") cpu="NA"
            if (wall=="") wall="NA"
            if (wall_sec=="") wall_sec="NA"
            if (mem_MB=="") mem_MB="NA"
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", file, tool, gene, user, sys, cpu, wall, wall_sec, mem_MB
        }
    ' "$f"
}

export -f parse_file

# Use GNU parallel if available, else fallback to serial
if command -v parallel >/dev/null 2>&1; then
    find "$logdir" -type f -name "*.time" | parallel -j$(nproc) parse_file {}
else
    for f in "$logdir"/*.time; do
        parse_file "$f"
    done
fi

