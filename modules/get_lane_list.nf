process get_lane_list {

    input:
    path complete_report

    output:
    path("${list_of_lanes}")

    script:
    list_of_lanes="lanes.txt"

    """
    # Create list of lanes that did not fail
    awk -F'\t' '\$0!~"FAIL" {print \$1}' ${complete_report} | sed '1d' > ${list_of_lanes}
    """
}
