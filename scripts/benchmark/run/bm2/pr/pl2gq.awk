#!/usr/bin/awk -f

# fix_het_gp is from cmdline

BEGIN { FS = "\t" }

{
    ra_same = $11;
    het_arr = $14;
    pre_gp = $7;

    n = split(pre_gp, a, ",");   # PL
    if (n != 3) exit 1;

    # get all GP
    for (i = 1; i <= n; i++) {
        s = 0;
        for (j = 1; j <= n; j++) {
            s += 10.0 ^ ((a[i] - a[j]) / 10.0);
        }
        gp[i] = 1.0 / s;
    }

    # get het GP
    if (fix_het_gp && het_arr && ! ra_same) {
        het_gp = 0;
    } else {
        het_gp = gp[2];
    }

    # get GQ
    m1 = a[2];
    m2 = a[1];
    if (a[1] < a[2]) {
        m1 = a[1];
        m2 = a[2];
    }
    for (i = 3; i <= n; i++) {
        if (a[i] < m1) {
            m2 = m1;
            m1 = a[i];
        } else if (a[i] < m2) {
            m2 = a[i];
        }
    }
    gq = m2 - m1;

    printf("%s\t%f\t%f,%f,%f\t%f\n", $0, gq, gp[1], gp[2], gp[3], het_gp);
}

#{
#    n = split($5, a, ",");
#    if (n < 3) exit 1;
#    if (a[1] < a[2]) {
#        m1 = a[1];
#        m2 = a[2];
#    } else { 
#        m1 = a[2];
#        m2 = a[1];
#    }
#    for (i = 3; i <= n; i++) {
#        if (a[i] < m1) {
#            m2 = m1;
#            m1 = a[i];
#        } else if (a[i] < m2) {
#            m2 = a[i];
#        }
#    }
#    s = 0;
#    for (i = 1; i <= n; i++) {
#        s += 10 ^ ((m1 - a[i]) / 10);
#    }
#    max_gp = 1 / s;
#    gq = m2 - m1;
#    printf("%s\t%f\t%f\n", $0, max_gp, gq);
#}
