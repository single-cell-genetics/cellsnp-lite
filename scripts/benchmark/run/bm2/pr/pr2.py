#Precision Recall Curve

import sys
import pylab as pl
import numpy
import matplotlib.pyplot as plt
import getopt

WIDTH = 4.5
HEIGHT = 3.5
DPI = 300

def PR_stat(state, scores, thresh):
    """
    @abstract       Non-standard Precision-Recall Curve
    @param state    If prediction is correct [vector <int>]
    @param scores   Scores [vector <float>]
    @param thresh   Marker threshold (single value), None if not use [float]
    @return         Stat results, a tuple of <precision> <recall> <threshold> <auc> if success,
                    None otherwise.
    """
    n = state.shape[0]
    if n <= 0 or scores.shape[0] != n:
        return None

    score_gap = numpy.array([thresh])
    if thresh is None:
        score_gap = numpy.unique(scores)
        if score_gap.shape[0] > 2000:
            idx = numpy.random.permutation(score_gap.shape[0])
            score_gap = score_gap[idx[:2000]]    
    thresholds = numpy.sort(score_gap)
    
    precision = numpy.zeros(thresholds.shape[0], dtype = "float")
    recall = numpy.zeros(thresholds.shape[0], dtype = "float")
    for i in range(thresholds.shape[0]):
        idx_r = numpy.where(scores >= thresholds[i])[0]
        nr = idx_r.shape[0]
        recall[i] = 1.0 * nr / n
        np = sum(state[idx_r])
        precision[i] = 1.0 * np / nr

    auc = None
    if thresh is None:
        _recall = numpy.append(1.0, recall)
        _recall = numpy.append(_recall, 0.0)
        _precision = numpy.append(0.0, precision)
        _precision = numpy.append(_precision, 1.0)
        auc = 0
        for i in range(_recall.shape[0] - 1):
            auc += (_recall[i] - _recall[i + 1]) * (_precision[i] + _precision[i + 1]) / 2.0

    return precision, recall, thresholds, auc

def PR_plot2(name1, infile1, name2, infile2, outfig, title, width, height, dpi):
    """
    @abstract      Plot Precision-Recall curve and calc AUC for genotyping results of one app
    @param name1   Name of app1 [str]
    @param infile1 Input tsv file of app1, two columns: <label> <prob> [str]
    @param name2   Name of app2 [str]
    @param infile2 Input tsv file of app2, two columns: <label> <prob> [str]
    @param outfig  Output figure [str]
    @param title   Title of figure [str]
    @param width   Width of figure in inch [float]
    @param height  Height of figure in inch [float]
    @param dpi     Dpi of figure [int]
    @return        Void
    """
    fig = pl.figure()

    for name, infile, color in zip((name1, name2), (infile1, infile2), ("red", "blue")):
        data = numpy.loadtxt(infile, delimiter = "\t")
        precision, recall, thresholds, auc = PR_stat(data[:, 0], data[:, 1], None)
        mpre, mrec, mth, mauc = PR_stat(data[:, 0], data[:, 1], 20)

        plt.plot(recall, precision, "-", color = color, linewidth = 1.0,
                label = "%s: AUC=%.3f" % (name, auc))
        plt.plot(mrec, mpre, marker = "o", markersize = 3, color = color)

    plt.plot(numpy.arange(0, 2), 1 - numpy.arange(0, 2), "k--", 
            linewidth=1.0, label="random: AUC=0.500")
    plt.legend(loc="best", fancybox=True, ncol=1)
    plt.xlabel("Recall")
    plt.ylabel("Precision")

    if title:
        pl.title(title)
    fig.set_size_inches(width, height)
    fig.savefig(outfig, dpi = dpi, bbox_inches = "tight")

def usage(fp = sys.stderr):
    msg = "\n"
    msg += "Usage: %s [options]\n" % sys.argv[0]
    msg += "\n"                                        \
           "Options:\n"                                   \
           "  --name1 STR      Name of app1: cellsnp-lite or freebayes.\n"   \
           "  --infile1 FILE   Input tsv file of app1\n"          \
           "  --name2 STR      Name of app2: cellsnp-lite or freebayes.\n"   \
           "  --infile2 FILE   Input tsv file of app2\n"          \
           "  --outfig FILE    Output figure\n"              \
           "  --title STR      Title of figure\n" 
    msg += "  --width FLOAT    Width of figure in inch [%f]\n" % WIDTH
    msg += "  --height FLOAT   Height of figure in inch [%f]\n" % HEIGHT 
    msg += "  --dpi INT        Dpi of figure [%d]\n" % DPI
    msg += "  -h, --help       Print this message\n"     \
           "\n"
    fp.write(msg)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["name1=", "infile1=", "name2=", 
                       "infile2=", "outfig=", "title=", "width=", "height=", "dpi=", "help"])
    except getopt.GetoptError as e:
        print(str(e))
        usage()
        sys.exit(2)
    name1 = infile1 = name2 = infile2 = outfig = title = width = height = dpi = None
    for opt, value in opts:
        if opt == "--name1": name1 = value
        elif opt == "--infile1": infile1 = value
        elif opt == "--name2": name2 = value
        elif opt == "--infile2": infile2 = value
        elif opt == "--outfig": outfig = value
        elif opt == "--title": title = value
        elif opt == "--width": width = float(value)
        elif opt == "--height": height = float(value)
        elif opt == "--dpi": dpi = int(value)
        elif opt in ("-h", "--help"): usage(); sys.exit(3)
        else: assert False, "unhandled option"
    if width is None: width = WIDTH
    if height is None: height = HEIGHT
    if dpi is None: dpi = DPI
    PR_plot2(name1, infile1, name2, infile2, outfig, title, width, height, dpi)

