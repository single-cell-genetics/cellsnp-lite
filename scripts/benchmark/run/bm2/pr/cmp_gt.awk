#!/usr/bin/awk -f

{
    ra_app = $3;
    gt_app = $5;
    tgt_app = $6;

    ra_arr = $8;
    gt_arr = $9;
    tgt_arr = $10;

    ra_same = 0;
    bi_same = 0;
    het_app = 0;
    het_arr = 0;

    n1 = split(ra_app, a, "|");    # ref_alt of app
    n2 = split(ra_arr, b, "|");    # ref_alt of array
    if (n1 != 2 || n2 != 2) { exit 1; }
    if ((a[1] == b[1] && a[2] == b[2]) || (a[1] == b[2] && a[2] == b[1])) {
        ra_same = 1;
    }

    n1 = split(tgt_app, a, "/");    # TGT of app
    n2 = split(tgt_arr, b, "/");    # TGT of array
    if (n1 != 2 || n2 != 2) { exit 1; }
    if ((a[1] == b[1] && a[2] == b[2]) || (a[1] == b[2] && a[2] == b[1])) {
        bi_same = 1;
    } 

    if (gt_app == "0/1" || gt_app == "1/0") { het_app = 1; }
    if (gt_arr == "0/1" || gt_arr == "1/0") { het_arr = 1; }

    printf("%s\t%d\t%d\t%d\t%d\n", $0, ra_same, bi_same, het_app, het_arr);
}
