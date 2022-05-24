import pandas as pd
from pyparsing import col


def get_statistics(bw_aln, bow_aln, bw_mis, bow_mis, bw_un, bow_un, bw_pre, bow_pre, bw_rec, bow_rec, bw_fs, bow_fs, path, number_of_sets):
    stats = [[], [], [], [], [], []]
    for i in range(number_of_sets):
        stats[0].extend([round(bw_aln[i], 2), round(bow_aln[i], 2)])
        stats[1].extend([round(bw_mis[i], 2), round(bow_mis[i], 2)])
        stats[2].extend([round(bw_un[i], 2), round(bow_un[i], 2)])
        stats[3].extend([round(bw_pre[i], 2), round(bow_pre[i], 2)])
        stats[4].extend([round(bw_rec[i], 2), round(bow_rec[i], 2)])
        stats[5].extend([round(bw_fs[i], 2), round(bow_fs[i], 2)])

    mux_header = ['Set ' + str(i+1) for i in range(number_of_sets)]
    mux = pd.MultiIndex.from_product([mux_header, ['BWA MEM','Bowtie']])
    stat_df = pd.DataFrame(stats, index=['Aligned', 'Misaligned', 'Unmapped', 'Precision', 'Recall', 'Fscore'], columns=mux)
    style_df = stat_df.style.set_table_styles([
                            {
                                "selector":"thead",
                                "props":"background-color:#F24BC0; color:white; border:3px;"
                            },
                            {
                                "selector":"th.row_heading",
                                "props":"background-color : #F24BC0; color:white; border:3px; text-align: left"}, 
                            {   "selector": "td", 
                                "props": "text-align:center; background-color: #F2F2F2"
                            },
                            
                            ])
    style_df.format(precision = 2)
    #style_df.to_html(path)
    stat_df.to_csv(path)
    
    


def set_errors_table():
        
    error_sets = {'Set 1' : [0, 0, 0], 'Set 2' : [0.0003, 0.0003, 0.0003], 'Set 3' : [0.03, 0, 0], 'Set 4' : [0, 0.03, 0.03], 'Set 5' : [0.03, 0.03, 0.03]}
    sets_df = pd.DataFrame(error_sets, index=['SNV', 'INS','DEL'])
    style_df = sets_df.style.set_table_styles([
                            {
                                "selector":"thead",
                                "props":"background-color:#F24BC0; color:white; border:3px;"
                            },
                            {
                                "selector":"th.row_heading",
                                "props":"background-color : #F24BC0; color:white; border:3px; text-align: left"}, 
                            {   "selector": "td", 
                                "props": "text-align:center; background-color: #F2F2F2"
                            },
                            
                            ])
    
    style_df.format(precision = 2, formatter = {'Set 2' : "{:.4f}"})
    style_df.to_html("./testing2/Error_set.html")

def genomes_info_table():
     genomes = { 'Human herpesvirus' : ["151 KB"], 'Clostridium tetani' : ['2.79 MB'], 'Mycobacterium tuberculosis' : ["4.25 MB"]}
     gen_df = pd.DataFrame(genomes, index=['Reference genome size'])
     style_df = gen_df.style.set_table_styles([
                            {
                                "selector":"thead",
                                "props":"background-color:#F24BC0; color:white; border:3px; font-style:italic;"
                            },
                            {
                                "selector":"th.row_heading",
                                "props":"background-color : #F24BC0; color:white; border:3px; text-align: left"}, 
                            {   "selector": "td", 
                                "props": "text-align:center; background-color: #F2F2F2"
                            },
                            
                            ])
     style_df.to_html("./testing2/Genomes_info.html")



#set_errors_table()
#genomes_info_table()