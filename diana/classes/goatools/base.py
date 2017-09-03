# Stolen from brentp:
# <https://github.com/brentp/toolshed/blob/master/toolshed/files.py>

import bz2
import gzip
import wget
import sys
import urllib
import os.path as op
import os


if sys.version_info[0] < 3:
    int_types = (int, long)
    urlopen = urllib.urlopen
else:
    int_types = (int,)
    basestring = str
    from urllib.request import urlopen


def nopen(f, mode="r"):
    r"""
    open a file that's gzipped or return stdin for '-'
    if f is a number, the result of nopen(sys.argv[f]) is returned.
    >>> nopen('-') == sys.stdin, nopen('-', 'w') == sys.stdout
    (True, True)
    >>> nopen(sys.argv[0])
    <...file...>
    # expands user and vars ($HOME)
    >>> nopen("~/.bashrc").name == nopen("$HOME/.bashrc").name
    True
    # an already open file.
    >>> nopen(open(sys.argv[0]))
    <...file...>
    >>> nopen(0)
    <...file...>
    Or provide nicer access to Popen.stdout
    >>> files = list(nopen("|ls"))
    >>> assert 'setup.py\n' in files or b'setup.py\n' in files, files
    """
    if isinstance(f, int_types):
        return nopen(sys.argv[f], mode)

    if not isinstance(f, basestring):
        return f
    if f.startswith("|"):
        # using shell explicitly makes things like process substitution work:
        # http://stackoverflow.com/questions/7407667/python-subprocess-subshells-and-redirection
        # use sys.stderr so we dont have to worry about checking it...
        p = Popen(f[1:], stdout=PIPE, stdin=PIPE,
                  stderr=sys.stderr if mode == "r" else PIPE,
                  shell=True, bufsize=-1, # use system default for buffering
                  preexec_fn=prefunc,
                  close_fds=False, executable=os.environ.get('SHELL'))
        if sys.version_info[0] > 2:
            import io
            p.stdout = io.TextIOWrapper(p.stdout)
            p.stdin = io.TextIOWrapper(p.stdin)
            if mode != "r":
                p.stderr = io.TextIOWrapper(p.stderr)

        if mode and mode[0] == "r":
            return process_iter(p, f[1:])
        return p

    if f.startswith(("http://", "https://", "ftp://")):
        fh = urlopen(f)
        if f.endswith(".gz"):
            return ungzipper(fh)
        if sys.version_info[0] < 3:
            return fh
        import io
        return io.TextIOWrapper(fh)
    f = op.expanduser(op.expandvars(f))
    if f.endswith((".gz", ".Z", ".z")):
        fh = gzip.open(f, mode)
        if sys.version_info[0] < 3:
            return fh
        import io
        return io.TextIOWrapper(fh)
    elif f.endswith((".bz", ".bz2", ".bzip2")):
        fh = bz2.BZ2File(f, mode)
        if sys.version_info[0] < 3:
            return fh
        import io
        return io.TextIOWrapper(fh)

    return {"r": sys.stdin, "w": sys.stdout}[mode[0]] if f == "-" \
         else open(f, mode)


def ungzipper(fh, blocksize=16384):
    """
    work-around to get streaming download of http://.../some.gz
    """
    import zlib
    uzip = zlib.decompressobj(16 + zlib.MAX_WBITS)
    data = uzip.decompress(fh.read(blocksize)).split("\n")

    while len(data[0]):
        # last chunk might not be a full line.
        save = data.pop()
        for line in data:
            yield line
        data = uzip.decompress(fh.read(blocksize)).split("\n")
        # first line is prepended with saved chunk from end of last set.
        data[0] = save + data[0]


def download_go_basic_obo(obo="go-basic.obo", prt=sys.stdout):
    """Download Ontologies, if necessary."""
    # Download: http://geneontology.org/ontology/go-basic.obo
    if not os.path.isfile(obo):
        slim = "subsets" if "slim" in obo else "."
        obo_remote = "http://geneontology.org/ontology/{SLIM}/{OBO}".format(SLIM=slim, OBO=obo)
        wget.download(obo_remote)
        if prt is not None:
            prt.write("\n  DOWNLOADED: {FILE}\n".format(FILE=obo))
    else:
        if prt is not None:
            prt.write("  EXISTS: {FILE}\n".format(FILE=obo))
    return obo

def download_ncbi_associations(gene2go="gene2go", prt=sys.stdout):
    """Download associations from NCBI, if necessary"""
    # Download: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    gz = "{GENE2GO}.gz".format(GENE2GO=gene2go)
    if not os.path.isfile(gene2go):
        wget.download("https://ftp.ncbi.nlm.nih.gov/gene/DATA/{GZ}".format(GZ=gz))
        assert gunzip(gz) == gene2go
        if prt is not None:
            prt.write("\n  DOWNLOADED: {FILE}\n".format(FILE=gene2go))
    else:
        if prt is not None:
            prt.write("  EXISTS: {FILE}\n".format(FILE=gene2go))
    return gene2go

def gunzip(gz, file_gunzip=None):
    """Unzip .gz file. Return filename of unzipped file."""
    if file_gunzip is None:
        file_gunzip = os.path.splitext(gz)[0]
    with gzip.open(gz, 'rb') as zstrm:
        with  open(file_gunzip, 'wb') as ostrm:
            ostrm.write(zstrm.read())
    os.remove(gz)
    return file_gunzip

def get_godag(fin_obo="go-basic.obo", prt=sys.stdout):
    """Return GODag object. Initialize, if necessary."""
    from goatools.obo_parser import GODag
    download_go_basic_obo(fin_obo, prt)
    return GODag(fin_obo)

def dnld_gaf(species_txt):
    """Download GAF file if necessary."""
    return dnld_gafs([species_txt])[0]

def dnld_gafs(species_list):
    """Download GAF files if necessary."""
    # Example GAF files:
    #   http://geneontology.org/gene-associations/gene_association.goa_human.gz
    #   http://geneontology.org/gene-associations/gene_association.mgi.gz
    #   http://geneontology.org/gene-associations/gene_association.fb.gz
    fin_gafs = []
    for species_txt in species_list: # e.g., goa_human mgi fb
        gaf_ftp = "http://geneontology.org/gene-associations/gene_association.{S}".format(
            S=species_txt)
        gaf_local = "gene_association.{SPECIES}".format(SPECIES=species_txt)
        if not os.path.isfile(gaf_local):
            wget.download("{GAF}.gz".format(GAF=gaf_ftp))
            gaf_gz = "{GAF_LOCAL}.gz".format(GAF_LOCAL=gaf_local)
            with gzip.open(gaf_gz, 'rb') as zstrm:
                with  open(gaf_local, 'wb') as ostrm:
                    ostrm.write(zstrm.read())
            os.remove(gaf_gz)
        fin_gafs.append(gaf_local)
    return fin_gafs

# Copyright (C) 2013-2017, B Pedersen, et al. All rights reserved."
