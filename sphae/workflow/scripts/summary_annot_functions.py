import pandas as pd


def create_summary(phold_func, pkl_func, tmp, summary_gbk):
    phold_df = pd.read_csv(phold_func, sep='\t', header=None)
    pkl_df   = pd.read_csv(pkl_func,   sep='\t', header=None)

    max_len = max(len(phold_df), len(pkl_df))
    phold_df = phold_df.reindex(range(max_len))
    pkl_df   = pkl_df.reindex(range(max_len))

    tmp_df = pd.concat([phold_df, pkl_df[2]], axis=1)
    tmp_df.to_csv(tmp, sep='\t', header=False, index=False)

    summary_df = pd.DataFrame({
        "contig name": phold_df[0],
        "protein ID":  phold_df[1],
        "phold":       phold_df[2],
        "phynteny":    pkl_df[2],
    })
    summary_df = summary_df.fillna("NA")
    summary_df.to_csv(summary_gbk, sep='\t', header=True, index=False)


create_summary(
    phold_func=snakemake.input.phold_func,
    pkl_func=snakemake.input.pkl_func,
    tmp=snakemake.params.tmp,
    summary_gbk=snakemake.output.summary_gbk,
)
