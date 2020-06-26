#!/usr/bin/env/ python
import sys
import os
import csv
import matplotlib.pyplot as plt
import seaborn as sns
import statistics

def open_csv_data(data_file):
    gRNA_counts = {}
    with open(data_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            # add line_count=0 +=1 to add header if needed
            gRNA_counts[row[0]] = int(row[1])
    return gRNA_counts

def library_description(gRNA_count):
    gRNA_counter = list(gRNA_count.values())
    minimum = min(gRNA_counter)
    maximum = max(gRNA_counter)
    mean = round(statistics.mean(gRNA_counter), 3)
    sd = round(statistics.stdev(gRNA_counter), 3)
    return minimum, maximum, mean, sd


def plot_gRNA(gRNA_count):
    sns.distplot(list(gRNA_count.values()), hist=True, kde=True,
                 # use hist=True to plot a histogram under the curve
                 # color='darkblue',
                 # hist_kws={'edgecolor': 'black'},
                 # kde_kws={'linewidth': 4}
                 )

def plot_violin_gRNA(gRNA_count):
    sns.violinplot(x=list(gRNA_count.values()))


try:
    filesToPlot = sys.argv[1:]
except IndexError:
    raise SystemExit(f"Usage: gRNA_dist_plot.py <files with gRNA counts separated by spaces>")

gRNA_count_data = []
gRNA_count_stats = []

for file in filesToPlot:
    gRNA_count_data.append(open_csv_data(file))


for sample in gRNA_count_data:
    PROVA = (filesToPlot[0],) + library_description(sample) #Iterate this zero
    gRNA_count_stats.append(PROVA)
    plot_gRNA(sample)
    #plot_violin_gRNA(sample) #it would be nice to make a violin plot, using pandas

print("minimum, maximum, mean, sd")
print(gRNA_count_stats)

legend_names = [os.path.splitext(name)[0] for name in filesToPlot] #removes .csv extension for legend
plt.legend(title='gRNA distribution', loc='upper right', labels=legend_names)
plt.savefig("gRNAdistribution.png")


# print("Min,Max,Mean,DesvEst Cas9P", library_description(gRNA_counts_Cas9P))
# print("Min,Max,Mean,DesvEst Cas9M", library_description(gRNA_counts_Cas9M))
# print("Min,Max,Mean,DesvEst CBE", library_description(gRNA_counts_CBE))
# print("Min,Max,Mean,DesvEst ABE", library_description(gRNA_counts_ABE))
# plot_gRNA(gRNA_counts_Cas9P)
# plot_gRNA(gRNA_counts_Cas9M)
# plot_gRNA(gRNA_counts_CBE)
# plot_gRNA(gRNA_counts_ABE)
# plt.legend(title='gRNA distribution', loc='upper right', labels=[plasmid_seq_Cas9P, plasmid_seq_Cas9M,
#                                                                  plasmid_seq_CBE, plasmid_seq_ABE])
# #plt.show()
# #plt.savefig('gRNA_distribution.png')
