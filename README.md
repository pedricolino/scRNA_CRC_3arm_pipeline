Originally, I just created a small SnakeMake pipeline to preprocess all samples individually according to and with Theis lab's best practices notebooks. 

Later, it evolved into a bigger pipeline to explore results when using different preprocessing parameters,  such as different methods to reduce the effect of cell cycle phase on the dimensional reductions and clustering, the optional use of a SoupX-corrected count matrix, and different quality control methods altogether (see scAutoQC). 

Additionally, I performed most analyses first interactively on a subset of the data to write and test scripts/notebooks, and then reran the same script/notebook on the full dataset via Snakemake.

In the end, I use jupyter-book to render those notebooks as a book/website to share.
