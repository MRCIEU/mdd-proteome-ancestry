# Using RAP

First attempt was to use swiss-army-knife. Use bgenix to extract a list of SNPs, then cat-bgen to concatenate. bgenix worked, but cat-bgen gives an error about the bgen file being wrong somehow. Plink is unable to read the bgen files, it says it has unexpected number of columns. But bgenix can convert to vcf and get info on them. Not clear what the issue is

Using ttyd with docker a bit frustrating because no obvious images to use for various tools.

