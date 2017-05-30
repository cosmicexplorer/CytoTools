auto-mem
========

**THIS FILE IS OLD AND USEFUL FOR MY (DANNY'S) NOTES, DON'T WORRY ABOUT IT**

# Overview

The `auto-mem.R` script performs MEM calculations between populations in FCS files. You can designate populations as 'references', and over every population, for each reference, the MEM of each channel shared with the reference is calculated and displayed along with other data. You can tell the script to remove or change the names of certain channels, and it cleans up the input files by default. It can be used on the command line or in an R interactive session.

You can run `auto-mem.R` from a command line or by loading it interactively and calling the `analyze()` function with an argument list (shown below).

# Requirements, Assumptions, and Behavior

**You need to install the following R packages before use:**

- `flowCore`
- `dplyr`
- `argparse`
- `stringr`
- `stringdist`
- `magrittr`

**The script will assume your fcs files each contain a single population.**

**The script will remove rows with NaNs, columns of all integers, and duplicate rows.**

# Example Usage

The following command is used to demonstrate the script's features and how to call it. This can be run from a command line.

``` bash
./auto-mem.R -v -f example-regexps-file.txt -f example-regexps-file-2.txt \
             -c example-canon-file.txt -c example-canon-file-2.txt \
             -r '1305007 D0 PB.fcs' -r '1305007 D8.fcs' \
             '1305007 D0 BM.fcs'
```

The above command does the following:

- `./auto-mem.R` = "*call* the auto-mem script"
- `-v` = "with *verbose* output"
- `-f example-regexps-file.txt -f example-regexps-file-2.txt` = "*exclude* columns with names matching any *patterns* from these file(s)"
- `-c example-canon-file.txt -c example-canon-file-2.txt` = "*canonicalize* column names as specified in these file(s)"
- `-r '1305007 D0 PB.fcs' -r '1305007 D8.fcs'` = "use these file(s) as *reference* populations in the analysis"
- `'1305007 D0 BM.fcs'` = "use these file(s) as *non-reference* populations in the analysis"

[Patterns](#patterns) and [canonicalization](#canonicalization) are described below. The above command corresponds to running the below script in an interactive R session (you *must* provide values for all of the keys listed below).

``` R
source("auto-mem.R")

print_text_output(analyze(list(
    column_matches_file = list("example-regexps-file.txt", "example-regexps-file-2.txt"),
    canonize_mapping_file = list("example-canon-file.txt", "example-canon-file-2.txt"),
    files = list("1305007 D0 BM.fcs"),
    ref = list("1305007 D0 PB.fcs", "1305007 D8.fcs"),
    verbose = T, no_strip_ints = T, use_perl_regexps = F,
    case_sensitive = F, channel_canon_input_file = NULL,
    channel_canon_output_file = NULL)))
```

# Input Files and Example Output

Let `patterns.txt` contain the following:

```
nd143
tm169
```

Let `canons.txt` contain the following:

```
(Er167)Dd:erddd
(Dy164)Dd:ddyydd
```

The output looks like this:

``` bash
$ ./auto-mem.R -v -f patterns.txt -c canons.txt \
             -r '1305007 D0 PB.fcs' -r '1305007 D8.fcs' \
             '1305007 D0 BM.fcs'
reading/processing fcs files...
columns removed by regexp (in file '1305007 D0 BM.fcs'): [(Nd143)Dd, (Tm169)Dd]
columns of all integers removed: [Time, Cell_length]
columns removed by regexp (in file '1305007 D0 PB.fcs'): [(Nd143)Dd, (Tm169)Dd]
columns of all integers removed: [Time, Cell_length]
columns removed by regexp (in file '1305007 D8.fcs'): [(Nd143)Dd, (Tm169)Dd]
columns of all integers removed: [Time, Cell_length]
(processing took 2.067000 seconds)
performing analysis...
(analysis took 0.506000 seconds)
1305007 D0 BM.fcs (abundance=24.67%, type=no_ref, file=1305007 D0 BM.fcs):
     (Ba137)Dd  (Ba138)Dd  (Pr141)Dd  (Nd142)Dd  (Nd143)Dd  (Nd144)Dd
25% -0.7102908 -0.5767109 -0.7155359 -0.6069692 -0.5890164 -0.5415837
50% -0.4185478 -0.1525823 -0.4293207 -0.2145401 -0.1804526 -0.0824318
75% -0.1261168  2.6756206 -0.1455588  1.7301667  2.8392251  3.9768801
... (more output)
```

# Glossary
## Patterns

Adding `-f COLUMN_MATCHES_FILE` to the command line reads the lines of `COLUMN_MATCHES_FILE` as patterns. By default, these are case-insensitive. Adding `-I` makes them case-sensitive and `-P` makes them PCRE regular expressions. The patterns are *partially matched*, so if they appear anywhere in the column name, the column will be removed.

For example, let `patterns.txt` be a file in the current directory with the following content.
```
tsne
visne
event
```

Running the command `./auto-mem.R -f patterns.txt ...` would remove these channel names from input fcs files (along with many more):

- `tSNE1`
- `data-visne-1`
- `event number`
- `Event#`

## Canonicalization

Adding `-c CANONIZE_MAPPING_FILE` to the command line reads the lines of `CANONIZE_MAPPING_FILE` as *canonicalization mappings*, which are read as `<channel name>:<canonical name>`. In the analysis, the *channel name* is transformed to *canonical name*.

### Restrictions

- **name matching is done case insensitively unless `-I` is given**
- **you cannot canonicalize a channel name with colons (`:`) in it!**

### Example Usage

For example, let `canons.txt` be a file in the current directory with the following content.
```
(Er167)Dd:erddd
(Dy164)Dd:ddyydd
```

Running the command `./auto-mem.R -c canons.txt ...` would transform *only* the following channel names:

- `(Er167)Dd` is changed to `erddd` (a random name I made up)
- `(Dy164)Dd` is changed to `ddyydd` (another random name I made up)

### When to Canonicalize

Canonicalization becomes useful when *multiple channel names refer to the exact same thing*. For example, say some of your fcs files contain a channel named `115-CD45`, others have `cd45 115`, and still others have `115 CD45`, and you know these refer to the same thing. A canonicalize file containing the following would let the script "understand" that these mean the same thing.

```
cd45 115:115-CD45
115 CD45:115-CD45
```

A file containg the above (if that file is given with `-c`) would let the script know that `115-CD45`, `cd45 115`, and `115 CD45` are the same channel, so it can compare these channels across populations. You could also use:

```
115-cd45:115_CD45_CANON
cd45 115:115_CD45_CANON
115 cd45:115_CD45_CANON
```

Because canonicalization is done *case insensitively* by default and the canonicalized name can be an existing name or not, this is the same as before **for the fcs files described above**.

### Canon Input/Output Files

The above used a toy example of canonicalization. The `-C CHANNEL_CANON_INPUT_FILE` and `-O CHANNEL_CANON_OUTPUT_FILE` arguments specify a file to read from and write to, respectively.

At start, the script will read lines of text from `CHANNEL_CANON_INPUT_FILE`. While processing fcs files, the script will report if any channel names used in the analysis are close to but not quite the same as a canonical name from `CHANNEL_CANON_INPUT_FILE` as well as other canonical names from the current set of fcs files.

If `CHANNEL_CANON_OUTPUT_FILE` is given, the resulting list of canonical names are written there, one per line, along with channel names that were NOT changed during canonicalization. This can be audited to see if any channels have not been correctly canonicalized.

You can specify the same file for `-C` and `-O`, but **be careful** when doing this unless you're sure you've selected the right files for each! The file will not be overwritten until after the canonization has been successfully applied (but before the analysis).

# Help and Usage

``` bash
$ ./auto-mem.R -h
usage: ./auto-mem.R [-h] [-v] [-f COLUMN_MATCHES_FILE] [-s] [-P] [-I] [-r REF]
                    [-c CANONIZE_MAPPING_FILE] [-C CHANNEL_CANON_INPUT_FILE]
                    [-O CHANNEL_CANON_OUTPUT_FILE]
                    FILE...

Process FCS files for MEM scores

positional arguments:
  FILE...               non-ref FCS files to be analyzed

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         print timing and debug output
  -f COLUMN_MATCHES_FILE, --column-matches-file COLUMN_MATCHES_FILE
                        remove columns with names matching patterns in file
  -s, --no-strip-ints   do NOT strip integer-only columns from input files
  -P, --use-perl-regexps
                        use PCRE regular expressions for patterns
  -I, --case-sensitive  make patterns case sensitive
  -r REF, --ref REF     use file as a reference population
  -c CANONIZE_MAPPING_FILE, --canonize-mapping-file CANONIZE_MAPPING_FILE
                        read channel canonicalization from this file
  -C CHANNEL_CANON_INPUT_FILE, --channel-canon-input-file CHANNEL_CANON_INPUT_FILE
                        file with one canonical channel name per line
  -O CHANNEL_CANON_OUTPUT_FILE, --channel-canon-output-file CHANNEL_CANON_OUTPUT_FILE
                        where to store all canonical channels from this analysis, in addition to
                        the -C argument (if supplied)

-c and -r can be specified multiple times to use multiple files. The script reads each line of
column match files as a pattern. By default, these are interpretened as fixed strings to be
partially matched, ignoring case. Adding -I stops ignoring case, and -P makes the script interpret
these as PCRE regular expressions. See http://www.regular-expressions.info/tutorial.html for a
tutorial, or https://regex101.com/ to test regular expressions.
```
