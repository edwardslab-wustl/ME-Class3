
import sys
import collections

import numpy as np
import sklearn.metrics
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib as mpl


class Sample:
    def __init__(self,sample_name, legend, predictions, expression, gene_ID):
        self.sample_name = sample_name
        #self.individual_sample_names = set(sample_name.split("_"))
        self.predictions = predictions
        self.expression = expression
        self.gene_ID = gene_ID    
        self.legend = legend

    def set_labels(self, labels):
        self.labels = labels

class Result:
    def __init__(self,sample_name,legend,labels,predictions,pred_scores):
        self.sample_name = sample_name
        self.legend = legend
        self.labels = labels
        self.predictions = predictions
        self.pred_scores = pred_scores

        self.roc_fpr = ()
        self.roc_tpr = ()
        self.roc_thresholds = ()

        self.acc_rejrate_accuracy = ()
        self.acc_rejrate_rejectRate = ()
        self.acc_rejrate_numGenes_P = ()
        self.acc_rejrate_numGenes_N = ()

#    def set_roc(self, fpr, tpr, thresholds):
#        self.roc_fpr = fpr
#        self.roc_tpr = tpr
#        self.roc_thresholds = thresholds

    def pull_roc(self):
        result_txt = "## " + self.legend + "," + self.sample_name + "\n"
        result_txt = result_txt + "#fpr,tpr\n"
        for i, fpr in enumerate(self.roc_fpr):
            tpr = self.roc_tpr[i]
            result_txt = result_txt + "\t".join((str(fpr),str(tpr))) + "\n"
        #print(result_txt)
        return result_txt

    def pull_acc_rejectrate(self):
        result_txt = "## " + self.legend + "," + self.sample_name + "\n"
        result_txt = result_txt + "#1-rejectRate,accuracy\n"
        for i, accuracy in enumerate(self.acc_rejrate_accuracy):
            rejectRate = self.acc_rejrate_rejectRate[i]
            result_txt = result_txt + "\t".join((str(rejectRate),str(accuracy))) + "\n"
        #print(result_txt)
        return result_txt

def set_matplotlib_params(args, version):
        # See https://matplotlib.org/users/customizing.html
        lineWidth = float(args.lineWidth)
        params = {}
        if version >= 2:
            params = {'legend.fontsize': 25,
                    'legend.frameon': False,
                    'legend.borderaxespad': 0.,
                    'legend.handlelength': 1,
                    'legend.handletextpad': 0.3,
                    'figure.figsize': (10, 10),
                    'lines.linewidth' : lineWidth,
                    'axes.linewidth' : 2,
                    'axes.labelsize': 30,
                    'axes.titlesize':'xx-large',
                    'axes.spines.top' : False,
                    'axes.spines.right' : False,
                    'xtick.direction' : 'out',
                    'ytick.direction' : 'out',
                    'xtick.top' : False,
                    'ytick.right' : False,
                    'xtick.major.size' : 15,
                    'ytick.major.size' : 15,
                    'xtick.major.width' : 2,
                    'ytick.major.width' : 2,
                    'xtick.labelsize': 25,
                    'ytick.labelsize': 25}
        else:
            warn1 = ("WARNING: Using matplot lib version < 2. Your are using " +
                "version ")
            warn2 = ("Everything should be okay, but no guarantees.\n")
            sys.stderr.write(warn1 + "%s.  " % mpl.__version__ + warn2)
            params = {'legend.fontsize': 25,
                    'legend.frameon': False,
                    'legend.borderaxespad': 0.,
                    'legend.handlelength': 1,
                    'legend.handletextpad': 0.3,
                    'figure.figsize': (10, 10),
                    'lines.linewidth' : lineWidth,
                    'axes.linewidth' : 2,
                    'axes.labelsize': 30,
                    'axes.titlesize':'xx-large',
                    'axes.spines.top' : False,
                    'axes.spines.right' : False,
                    'xtick.direction' : 'in',
                    'ytick.direction' : 'in',
                    #'xtick.top' : False,
                    #'ytick.right' : False,
                    'xtick.major.size' : 15,
                    'ytick.major.size' : 15,
                    'xtick.major.width' : 2,
                    'ytick.major.width' : 2,
                    'xtick.labelsize': 25,
                    'ytick.labelsize': 25}
        pylab.rcParams.update(params)
        return 0

def print_acc_rejectrate (results, type, args):
    outFile = args.outFileBase + "." + type + ".acc_rejectrate.txt"
    sys.stderr.write("writing Accuracy and RejectRate data to: %s\n" % outFile)
    fh = open(outFile, 'w')
    for result in results:
        fh.write(result.pull_acc_rejectrate())
    fh.close()
    return 0

def plot_acc_rejectrate (results, outFile, args,version):
    sys.stderr.write("plotting Accuracy and RejectRate data to: %s\n" % outFile)
    fig = plt.figure()
    for result in results:
        plt.plot(result.acc_rejrate_rejectRate, result.acc_rejrate_accuracy, label=result.legend)
    plt.xlabel("1 - Reject Rate")
    plt.ylabel("Accuracy")
    axes=plt.gca()
    axes.set_xlim([0.,1.])
    axes.set_ylim([0.,1.])
    if version < 2:
        axes.yaxis.set_ticks_position('left')
        axes.xaxis.set_ticks_position('bottom')
    legendLoc = args.legendLoc
    if legendLoc == "right":
        lgd=plt.legend(bbox_to_anchor=(1.1, 1),loc=2)
    else:
        lgd=plt.legend(loc="best")
    #plt.tight_layout()
    plt.savefig(outFile, bbox_extra_artists=(lgd,), bbox_inches='tight')

def print_roc (results, args):
    outFile = args.outFileBase + ".roc.txt"
    sys.stderr.write("writing ROC data to: %s\n" % outFile)
    roc_fh = open(outFile, 'w')
    for result in results:
        roc_fh.write(result.pull_roc())
    roc_fh.close()
    return 0

def plot_roc (results, outFile, args, version):
    sys.stderr.write("plotting ROC data to: %s\n" % outFile)
    fig = plt.figure()
    for result in results:
        plt.plot(result.roc_fpr, result.roc_tpr, label=result.legend)
    plt.xlabel("FPR")
    plt.ylabel("TPR")
    axes=plt.gca()
    axes.set_xlim([0.,1.])
    axes.set_ylim([0.,1.])
    axes.spines['right'].set_visible(True)
    axes.spines['top'].set_visible(True)
    if version < 2:
        axes.yaxis.set_ticks_position('left')
        axes.xaxis.set_ticks_position('bottom')
    legendLoc = args.legendLoc
    if legendLoc == "right":
        lgd=plt.legend(bbox_to_anchor=(1.1, 1),loc=2)
    else:
        lgd=plt.legend(loc="best")
    plt.savefig(outFile, bbox_extra_artists=(lgd,), bbox_inches='tight')
    return 0

def calculate_acc_rejectrate (results, type, args):
    sys.stderr.write("calculating accuracy and reject rate, class: %s\n" % type)
    steps = args.acc_rejectrate_steps
    for result in results:
        if args.verbose:
            sys.stderr.write("\ton sample %s\n" % (result.sample_name))
        accuracy_list = list()
        rejectRate_list = list()
        numGenes_P_list = list()
        numGenes_N_list = list()
        for threshold in np.linspace(0,0.5,steps):
            #sys.stderr.write("\t\tthresh %s\n" % (str(threshold)))
            TP_P = 0
            TP_N = 0
            numGenes_P = 0
            numGenes_N = 0
            totalGenes = 0
            totalGenes_P = 0
            totalGenes_N = 0
            for i,prediction in enumerate(result.pred_scores):
                totalGenes += 1
                if result.labels[i] < 0:
                    totalGenes_N += 1
                else:
                    totalGenes_P += 1
                if prediction >= 1-threshold:
                    numGenes_P += 1
                    if result.labels[i] == result.predictions[i]:
                        TP_P += 1
                elif prediction <= threshold:
                    numGenes_N += 1
                    if result.labels[i] == result.predictions[i]:
                        TP_N += 1
            accuracy = 0
            rejectRate = 0
            if type == 'all':
                if (numGenes_P + numGenes_N) > 0:
                    accuracy = float(TP_P + TP_N) / float(numGenes_P + numGenes_N)
                if totalGenes > 0:
                    rejectRate = float(numGenes_P + numGenes_N) / float(totalGenes)
            elif type == "pos":
                if (numGenes_P) > 0:
                    accuracy = float(TP_P) / float(numGenes_P)
                if totalGenes_P > 0:
                    rejectRate = float(numGenes_P) / float(totalGenes_P)
            elif type == "neg":
                if (numGenes_N) > 0:
                    accuracy = float(TP_N) / float(numGenes_N)
                if totalGenes_N > 0:
                    rejectRate = float(numGenes_N) / float(totalGenes_N)
            else:
                sys.stderr.write("\tinvalid type in calculate_acc_rejectrate: %s\n" 
                        % (type))
                exit()
            accuracy_list.append(accuracy)
            rejectRate_list.append(rejectRate)
            numGenes_P_list.append(numGenes_P)
            numGenes_N_list.append(numGenes_N)
        result.acc_rejrate_accuracy = accuracy_list
        result.acc_rejrate_rejectRate = rejectRate_list
        result.acc_rejrate_numGenes_P = numGenes_P_list
        result.acc_rejrate_numGenes_N = numGenes_N_list
    return 0


def calculate_roc (results, args):
    sys.stderr.write("calculating roc\n")
    for result in results:
        if args.verbose:
            sys.stderr.write("\ton sample %s\n" % (result.sample_name))
        fpr, tpr, thresholds = sklearn.metrics.roc_curve(result.labels,
                result.pred_scores)
        result.roc_fpr = fpr
        result.roc_tpr = tpr
        result.roc_thresholds = thresholds

    return 0

def calculate_labels (samples, args):
    sys.stderr.write("calculating labels\n")
    for sample in samples:
        if args.verbose:
            sys.stderr.write("\ton sample %s\n" % (sample.sample_name))
        labels = list()
        for expr in sample.expression:
            if expr > 0:
                labels.append(1)
            else:
                labels.append(-1)
        sample.set_labels(labels)
        #sys.stderr.write("\texample labels %s\n" % (str(labels[0])))
    return 0


def calculate_results(samples, args):
    results = list()
    sys.stderr.write("calculating results\n")
    for sample in samples:
        if args.verbose:
            sys.stderr.write("\ton sample %s\n" % (sample.sample_name))
        predictions = list()
        for pred in sample.predictions:
            if pred >= 0.5:
                predictions.append(1)
            else:
                predictions.append(-1)
        results.append( Result(sample.sample_name,sample.legend,sample.labels,
            predictions,sample.predictions) )
    return results
    
def print_label_counter(Y):
    counter = collections.Counter(Y)
    for k in sorted(counter.keys()):
        sys.stderr.write("Class: %s, #: %d\n" % (str(k),counter[k]))
    sys.stderr.write("\n")
    return

def get_function_name(func):
    return str(func).split('(')[0]