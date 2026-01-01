README
================

## Course Information

BIO709 Biostatistics

Instructor - Akira Terui

## Suggested Repository Structure

Github repository should be organized so you can locate files easily. My
suggest

- `/code` : may contain script files (`.R`).

- `/data_raw` : may contain raw data files. This data should not be
  edited and

- `/data_fmt` : may contain formatted or derived data (e.g., after data
  cleaning etc.).

I do NOT recommend creating sub-sub-directories under any of the above
sub-directories. It will create issues when you call files in those deep
directories, as you must specify `sub-dir/sub-sub-dir/sub-sub-sub...`. I
would rather set file naming rules within each sub-directory so that
files are sorted by file type (see **File Name** section).

## File name

- Don’t

  - **Don’t use white-space, comma, and period.** Never, ever. Use
    underscore instead (`_`)

  - Avoid **capital** in file names.

- Do

  - Use **pre-fix** to group files automatically. For example, you may
    use `figure_xxx.R` for files generating figures, `data_xxx.csv` for
    data files, etc. This will help you organize files because (1) they
    are automatically grouped together within a sub-directory, and (2)
    you can call specific types of files using `list.files()` in R.

  - Be mindful about the hierarchy when naming files. For example, the
    prefix may begin from the largest grouping category (e.g., figure,
    table, analysis, format, etc), then specify specific tasks within
    each file group.
