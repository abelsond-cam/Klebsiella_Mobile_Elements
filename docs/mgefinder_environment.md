# MGEfinder environment

## Check bwa and bowtie2 are installed

```bash
micromamba run -n mgefinder_env bwa
micromamba run -n mgefinder_env bowtie2-build --version
```
If either command is not found, reinstall (they can disappear after pip installs in the same env):
```bash
micromamba install -n mgefinder_env -c bioconda -c conda-forge bwa bowtie2
```
Re-check after any pip changes to the env.

## 1. numpy/scipy binary incompatibility

If mgefinder fails with:
```text
ValueError: numpy.dtype size changed, may indicate binary incompatibility. Expected 96 from C header, got 88 from PyObject
```
do:

```bash
micromamba run -n mgefinder_env pip install --force-reinstall numpy scipy
micromamba run -n mgefinder_env pip install mgefinder --no-deps
```

Pip may report "dependency conflicts"; that is expected with `--no-deps` and the install still succeeds.

## 2. Bio.Alphabet removed (Biopython too new)

If mgefinder fails with:
```text
ImportError: Bio.Alphabet has been removed from Biopython.
```
your env has Biopython ≥1.78; MGEfinder 1.0.6 needs an older Biopython that still has `Bio.Alphabet`. Pin Biopython to a version before 1.78:

```bash
micromamba run -n mgefinder_env pip install "biopython>=1.75,<1.78"
```

## 3. AttributeError: 'str' object has no attribute 'decode' (mgefinder dependency check)

If mgefinder fails during "CHECKING DEPENDENCIES" with:
```text
AttributeError: 'str' object has no attribute 'decode'
```
the package’s dependency checker is written for Python 2; under Python 3, `shell()` can return a str, so `.decode('utf-8')` fails. Patch the installed file:

1. Find it:
   ```bash
   micromamba run -n mgefinder_env python -c "import mgefinder.dependencies as d; print(d.__file__)"
   ```
2. Open that file and find the line (around line 12) that looks like:
   `output = shell(cmd.format(tool=self.tool), read=True).decode('utf-8').strip()`
3. Replace it with:
   ```python
   raw = shell(cmd.format(tool=self.tool), read=True)
   output = (raw.decode('utf-8') if isinstance(raw, bytes) else raw).strip()
   ```
4. Save and re-run the pipeline.

## 4. RuntimeError: generator raised StopIteration (formatbam)

If mgefinder formatbam fails with:
```text
File ".../mgefinder/formatbam.py", line 16, in read_sam_pairs
    yield(next(samfile), next(samfile))
StopIteration
RuntimeError: generator raised StopIteration
```
the generator raises `StopIteration` when the SAM runs out of read pairs; in Python 3.7+ that becomes `RuntimeError`. Patch the installed **formatbam.py**:

1. Find it:
   ```bash
   micromamba run -n mgefinder_env python -c "import mgefinder.formatbam as f; print(f.__file__)"
   ```
2. Find the function `read_sam_pairs` (around lines 16–18). Replace:
   ```python
   def read_sam_pairs(samfile):
    while True:
        yield(next(samfile), next(samfile))
   ```
   with:
   ```python
   def read_sam_pairs(samfile):
       while True:
           try:
               r1 = next(samfile)
               r2 = next(samfile)
           except StopIteration:
               return
           yield (r1, r2)
   ```
3. Save and re-run the pipeline.

If the error persists, the SAM may have no (or an odd number of) alignments left after removing secondaries; check the BWA input and alignment quality.

## 5. KeyError in pair step: columns not in index

If mgefinder **pair** fails with:
```text
KeyError: "['pair_id', 'spanning_count', 'direct_repeat_reference'] not in index"
```
the code selects the final output header at several steps, but `pair_id`, `spanning_count`, and `direct_repeat_reference` are added in *later* steps, so they are missing and cause KeyError. Patch the installed **pair.py** to select only columns that exist at each step:

1. Find it:
   ```bash
   micromamba run -n mgefinder_env python -c "import mgefinder.pair as p; print(p.__file__)"
   ```
2. Apply these three edits (use “select only columns that exist” in each case):

   - **assign_pairs** (around line 245): change
     `assigned_pairs = sorted_pairs.query("keep_pair == True").loc[:, self.get_header_list()].sort_values(...)`
     to:
     ```python
     _cols = sorted_pairs.columns.intersection(self.get_header_list())
     assigned_pairs = sorted_pairs.query("keep_pair == True").loc[:, _cols].sort_values(['contig', 'pos_5p', 'pos_3p'])
     ```
   - **count_insertion_spanning_reads** (around line 269): change
     `assigned_pairs = assigned_pairs.loc[:, self.get_header_list()]`
     to:
     ```python
     assigned_pairs = assigned_pairs.loc[:, assigned_pairs.columns.intersection(self.get_header_list())]
     ```
   - **filter_junction_spanning** (around line 254): change
     `pairs = pairs[self.get_header_list()]`
     to:
     ```python
     pairs = pairs[pairs.columns.intersection(self.get_header_list())]
     ```
   - **get_direct_repeats** (around line 280): change
     `flank_pairs = flank_pairs.drop(['direct_repeat_reference'], axis=1).merge(positions, how='left')`
     to:
     ```python
     if 'direct_repeat_reference' in flank_pairs.columns:
         flank_pairs = flank_pairs.drop(['direct_repeat_reference'], axis=1)
     flank_pairs = flank_pairs.merge(positions, how='left')
     ```

4. **inferseq_assembly** (pandas list assignment fix):
   In `mgefinder/inferseqassembly.py` (around line 45), change:
   ```python
   inferred_sequences_with_context[['loc']] = loc
   ```
   to:
   ```python
   inferred_sequences_with_context['loc'] = loc
   ```

5. Save and re-run the pipeline.

Optionally, also install pandas 1.x via micromamba so the rest of the pipeline behaves as expected:
```bash
micromamba install -n mgefinder_env -c conda-forge "pandas>=1.0,<2"
```
