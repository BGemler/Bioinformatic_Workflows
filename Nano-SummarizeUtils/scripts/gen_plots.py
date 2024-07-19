import csv
import seaborn as sns
import pandas as pd

# set matplotlib cache to writeable path
import os
os.environ['MPLCONFIGDIR'] = "/tmp/"
import matplotlib.pyplot as plt


def set_blank_to_zero(count):
    if count == "":
        count = 0
    else:
        count = int(float(count))

    return count

def gen_count_boxplot(count_loc, plt_out_loc):
    """
    """
    with open(count_loc, "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        sample_count_data_dict = {}
        sample_count_data_dict["sample_ids"] = []
        sample_count_data_dict["sample_types"] = []
        sample_count_data_dict["count_type"] = []
        sample_count_data_dict["count"] = []

        for row in reader:
            sample_id, raw, trim, hostr, k2_prot = row[:5]

            # correct counts to ints, set blanks to 0
            raw, trim, hostr, k2_prot = set_blank_to_zero(raw), set_blank_to_zero(trim), \
                                    set_blank_to_zero(hostr), set_blank_to_zero(k2_prot)


            # determine whether sample_ID is pos control, neg control, or experimental
            if "libraryblank" in sample_id.lower():
                sample_type = "Library_Blank"
            elif "ntc" in sample_id.lower():
                sample_type = "NTC"
            elif "zymostd" in sample_id.lower():
                sample_type = "Pos_Control"
            else:
                sample_type = "Experimental"
            
            sample_count_data_dict["sample_ids"].append(sample_id)
            sample_count_data_dict["sample_types"].append(sample_type)
            sample_count_data_dict["count_type"].append("raw")
            sample_count_data_dict["count"].append(raw)

            sample_count_data_dict["sample_ids"].append(sample_id)
            sample_count_data_dict["sample_types"].append(sample_type)
            sample_count_data_dict["count_type"].append("trim")
            sample_count_data_dict["count"].append(trim)

            sample_count_data_dict["sample_ids"].append(sample_id)
            sample_count_data_dict["sample_types"].append(sample_type)
            sample_count_data_dict["count_type"].append("hostr")
            sample_count_data_dict["count"].append(hostr)

            sample_count_data_dict["sample_ids"].append(sample_id)
            sample_count_data_dict["sample_types"].append(sample_type)
            sample_count_data_dict["count_type"].append("k2_prot")
            sample_count_data_dict["count"].append(k2_prot)
    f.close()

    # Build seaborn paired boxplot of data
    sns.set_theme(style="ticks")

    df_for_boxplot = pd.DataFrame(data = sample_count_data_dict)

    box = sns.boxplot(x="count_type", y="count",
                hue="sample_types", palette=["m", "g"],
                data=df_for_boxplot)
    box.set_yscale("log")

    # horizontal line at 1 million reads
    box.axhline(y = 1e6, xmin = 0, xmax = 1, \
                    color = "black", linestyle = "dashed")

    plt.savefig(plt_out_loc)
    plt.close()

    return


def gen_domain_count_barplot(sample_domain_taxid_counts, plt_out_loc):
    """
    """
    # First, get the order of sample_ids correct
    type_sampleid_dict = {}
    for sample_id in sample_domain_taxid_counts:
        # determine whether sample_ID is pos control, neg control, or experimental
        if "libraryblank" in sample_id.lower():
            sample_type = "Library_Blank"
        elif "ntc" in sample_id.lower():
            sample_type = "NTC"
        elif "zymostd" in sample_id.lower():
            sample_type = "Pos_Control"
        else:
            sample_type = "Experimental"
        
        if sample_type not in type_sampleid_dict:
            type_sampleid_dict[sample_type] = []
        type_sampleid_dict[sample_type].append(sample_id)

    sample_types = []
    plot_Sample_ids = []
    domain_count_dict = {}
    for sample_type in type_sampleid_dict:
        sample_ids = type_sampleid_dict[sample_type]

        for sample_id in sample_ids:
            sample_id_domain_counts = sample_domain_taxid_counts[sample_id]
            sample_types.append(sample_type)
            plot_Sample_ids.append(sample_id)

            for domain in sample_id_domain_counts:
                count = sample_id_domain_counts[domain]
                if domain not in domain_count_dict:
                    domain_count_dict[domain] = 0
                domain_count_dict[domain] += count
        
    # sort domains in order from most to least abundant
    domain_counts = []
    for domain in domain_count_dict:
        domain_counts.append([domain, domain_count_dict[domain]])
    domain_counts = sorted(domain_counts, key = lambda x:x[1], reverse = True)

    # intialize domain-sampleid-ordered-dict
    ordered_domains = []
    domain_Sampleid_ordered_counts = {}
    for domain, _ in domain_counts:
        ordered_domains.append(domain)
        domain_Sampleid_ordered_counts[domain] = []

    # Now, create domain-specific list
    # Create them with fractions, not integers, for normalizing!
    for sample_id in sample_domain_taxid_counts:
        sample_id_domain_counts = sample_domain_taxid_counts[sample_id]
        sum_count = sum(sample_id_domain_counts.values())

        for domain in ordered_domains:
            if domain in sample_id_domain_counts:
                count = sample_id_domain_counts[domain] / sum_count
            else:
                count = 0.0
            
            domain_Sampleid_ordered_counts[domain].append(count)
    
    # Create dataframe
    df_for_barplot = pd.DataFrame(data = domain_Sampleid_ordered_counts, index = plot_Sample_ids)

    # Build seaborn paired boxplot of data
    sns.set_theme(style="ticks")

    # .groupby(df_for_barplot.index.name)
    df_for_barplot.plot(kind = 'bar', stacked = True, width = 1.0)
    
    plt.xlabel("Samples")
    plt.ylabel("Domain Composition")
    plt.ylim([0, 1])
    plt.legend(bbox_to_anchor=(1.1, 1.0))

    plt.xticks(rotation=90)

    plt.savefig(plt_out_loc, bbox_inches="tight")
    plt.close()

    return
