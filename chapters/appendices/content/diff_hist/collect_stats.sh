#!/usr/bin/env bash

if [ "$#" -ne 1 ]; then
	echo "Please provide a directory as an argument."
	exit -1
fi

parent_dir=$1

for stats_path in $parent_dir/snap061_plots/{stats_*,perc_diff_*}; do
	stats_file=$(basename "$stats_path")
	out_file=${stats_file/_\(/_allsnaps_\(}
	echo "Merging stats for $stats_file..."

	for snap_dir in $parent_dir/snap*_plots; do
		if [ -e $snap_dir/$stats_file ]; then
			snap_num=$(basename "$snap_dir")
			echo -n "${snap_num:5:2}  "
			cat $snap_dir/$stats_file | cut -d' ' -f 2-
		fi
	done | column -t > $parent_dir/$out_file

	echo "Stats written to $out_file..."

done

