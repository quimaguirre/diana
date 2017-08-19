import sys
import argparse
import os
import biana
from biana import *


def main():

    options = parse_user_arguments()
    translate(options)


def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Translate a bPPI network into your selected codes",
        epilog      = "@oliva's lab 2014")
    parser.add_argument('-i','--network_input_file',dest='input',action = 'store',default='edges.txt',
                        help = 'Input file (default is edges.txt)')
    parser.add_argument('-n','--nodes_input_file',dest='input_nodes',action = 'store',default='nodes.txt',
                        help = 'Input file with nodes (default is nodes.txt)')
    parser.add_argument('-trans','--translation_of_nodes_file',dest='translation_file',action = 'store',default=None,
                        help = 'File with the translation of codes from BIANA to the selected type for all nodes')
    parser.add_argument('-iformat','--input_format',dest='iformat',action = 'store',default='guild',
                        help = 'Format of input files of edges:\tguild (default) or \tnetscore')
    parser.add_argument('-oformat','--output_format',dest='oformat',action = 'store',default='guild',
                        help = 'Format of output files of edges:\tguild (default) or \tnetscore')
    parser.add_argument('-all','--use_all_nodes',dest='alln',action = 'store_true',
                        help = 'Flag to use all nodes even without translation of codes')
    parser.add_argument('-x','--exclude_geneID',dest='xgid',action = 'store_true',default=True,
                        help = 'Exclude geneID codes as input (if you use BIANA codes as input you cannot mix them with geneID)')
    parser.add_argument('-xf','--force_non_BIANA',dest='fgid',action = 'store_true',default=False,
                        help = 'Force codes different than BIANA entities as input (default is false)')
    parser.add_argument('-taxid','--TaxID',dest='taxid',action = 'store',default=None,
                        help = 'Tax ID (i.e. human=9606 ) only applied if GeneSymbol codes')
    parser.add_argument('-ttype','--translation_type',dest='ttype',action = 'store',default='uniprotentry',
                        help = '''Type of identifier for the output translation of codes (default is uniprotentry )
                        Using "proteinsequence" provides with the longest sequence of all codes''')
    parser.add_argument('-b','--biana_version',dest='biana_version',action = 'store',default='BIANA_MARCH_2013',
                        help = 'BIANA version MySQL')
    parser.add_argument('-score','--Non-seed_score',dest='score',action = 'store',default='0.0',type=float,
                        help = 'Score of non-seed nodes (default is 0.0)')
    parser.add_argument('-o','--output_file',dest='output',action = 'store', default=sys.stdout,
                        help = 'Output file with edges(default is standard output)')
    parser.add_argument('-v','--verbose',dest='show',action = 'store_true',
                        help = 'Verbose execution')

    options=parser.parse_args()

    return options
def printverbose(f,flag,message):
 """Define verbose to print in file 'f', if 'flag'=True, message given as 'message'"""
 if flag: f.write("%s"%(message))


def printfasta (out,name,seq):
   out.write(">%s\n"%(name))
   n=0;
   w=""
   for s in seq:
    w += s
    n+=1
    if n == 80:
     out.write("%s\n"%(w))
     n=0
     w=""
   if n>0: out.write("%s\n"%(w))

def longestseq(uE):
   seq=""
   for s in uE.get_attribute(attribute_identifier='proteinsequence'):
    if len(s.value.get_sequence().upper())>=len(seq):
     seq=s.value.get_sequence().upper()
   return seq


def fileExist (file):               #Checks if a file exists AND is a file
 if file is not None: return os.path.exists(file) and os.path.isfile(file)
 else: return False

def translate(options):
    log=sys.stderr
    if not fileExist(options.input):
     print("File with input network is missing\n")
     sys.exit(10)

    new={}
    if fileExist(options.translation_file):
     ft=open(options.translation_file,"r")
     for line in ft:
      word=line.strip().split("\t")
      node=word[0]
      new.setdefault(node,set())
      for trans in word[1].split("'"):
       if trans != "" and trans != ",":
         name="_".join([str(x) for x in trans.split()])
         new[node].add(name)
     ft.close()



    if options.xgid:
     id_type=['accessionnumber','uniprotentry','genesymbol','pdb','refseq','uniprotaccession',
            'proteinsequence','unigene']
    elif options.fgid:
     id_type=['geneid','accessionnumber','uniprotentry','genesymbol','pdb','refseq','uniprotaccession',
            'proteinsequence','unigene']
    else:
     id_type=['accessionnumber','uniprotentry','genesymbol','pdb','refseq','uniprotaccession',
            'geneid','proteinsequence','unigene']

    fd=open(options.input,"r")
    file_network=options.output+".edges"
    file_nodes=options.output+".nodes"
    if not options.output == sys.stdout:
        out_network=open(file_network,"w")
        out_nodes=open(file_nodes,"w")
    else:
        out_network=sys.stdout
        out_nodes=sys.stdout
    code=set()
    translation={}
    for line in fd:
      if options.iformat == 'netscore' :
       word=line.split()
       p=word[0].upper()
       q=word[1].upper()
       score=float(word[2])
       code.add(p)
       code.add(q)
      else:
       word=line.split()
       p=word[0].upper()
       score=float(word[1])
       q=word[2].upper()
       code.add(p)
       code.add(q)
    fd.close()
    if len(code) <= 0 :
     raise IOError("Input File is empty")
    session = create_new_session(sessionID="biana_session",
                               dbname=options.biana_version,
                               dbhost="sbi.upf.edu",
                               dbuser="biana_user",
                               dbpassword="biana_password",
                               unification_protocol="uniprot_geneID_seqtax")
    lost_and_found=set()
    if options.fgid:
     for name in code:
      printverbose(log,options.show,("Code %s is not integrated by BIANA - GID(%s) -\n"%(name,options.fgid)))
      lost_and_found.add(name)
    else:
     for name in code:
      try:
       uE = session.get_user_entity(name)
       seq=longestseq(uE)
       printverbose(log,options.show,("Code %s found by BIANA code %s\n"%(name,repr(uE))))
      except Exception as e:
       print("Error: %s when searching %s\n" %(e,name))
      if len(seq)>1 :
       translation[name]=uE
      else:
       printverbose(log,options.show,("Code %s is not integrated by BIANA - GID(%s) -\n"%(name,options.fgid)))
       lost_and_found.add(name)

    for name in lost_and_found:
     found=False
     for typ in id_type:
       printverbose(log,options.show,("Trying id_type %s\n"%(typ)))
       if not found:
        if (typ == 'genesymbol' and options.taxid is not None):
         proteome = session.create_new_user_entity_set(
                   identifier_description_list =[name],
                   attribute_restriction_list=[("taxid",options.taxid)],
                   id_type=typ,
                   new_user_entity_set_id="proteome")
        else:
         proteome = session.create_new_user_entity_set(
                   identifier_description_list =[name],
                   id_type=typ,
                   new_user_entity_set_id="proteome")
        try:
         for i in proteome.get_user_entity_ids():
          uE = session.get_user_entity(i)
          printverbose(log,options.show,("BIANA id code %s for %s\n"%(i,name)))
          seqi=longestseq(uE)
          if len(seqi)>1:
           translation[name]=uE
           seq=seqi
           found=True
          if len(seqi)>len(seq):
           translation[name]=uE
           found=True
           seq=seqi
        except:
         found=False
     if not found:
        printverbose(log,options.show,("Code %s not found in BIANA Database\n"%(name)))
        new.setdefault(name,set()).add(name)

    for name in translation:
      uE=translation[name]
      for current_id in uE.get_attribute(attribute_identifier=options.ttype):
        new.setdefault(name,set()).add(current_id.value.upper())
    nodes=set()
    fd=open(options.input,"r")
    for line in fd:
      if options.iformat == 'netscore' :
       word=line.split()
       p=word[0].upper()
       q=word[1].upper()
       score=float(word[2])
      else:
       word=line.split()
       p=word[0].upper()
       score=float(word[1])
       q=word[2].upper()
      info="".join([str(word[i]) for i in range(3,len(word))])
      try:
       if len(new[p])>0:
        for a in new[p]:
         if  len(new[q])>0:
          for b in new[q]:
           nodes.add(a)
           nodes.add(b)
           if options.oformat == 'netscore' :
             out_network.write("\t{0}\t{1}\t{2:f}\t\t{3:s}\n".format(a,b,score,info))
           else:
             out_network.write("{0} {1:f} {2}\t\t{3:s}\n".format(a,score,b,info))
         elif options.alln:
           b=q
           nodes.add(a)
           nodes.add(b)
           if options.oformat == 'netscore' :
             out_network.write("\t{0}\t{1}\t{2:f}\t\t{3:s}\n".format(a,b,score,info))
           else:
             out_network.write("{0} {1:f} {2}\t\t{3:s}\n".format(a,score,b,info))
       elif options.alln:
         a=q
         if  len(new[q])>0:
          for b in new[q]:
           nodes.add(a)
           nodes.add(b)
           if options.oformat == 'netscore' :
             out_network.write("\t{0}\t{1}\t{2:f}\t\t{3:s}\n".format(a,b,score,info))
           else:
             out_network.write("{0} {1:f} {2}\t\t{3:s}\n".format(a,score,b,info))
         elif options.alln:
           b=q
           nodes.add(a)
           nodes.add(b)
           if options.oformat == 'netscore' :
             out_network.write("\t{0}\t{1}\t{2:f}\t\t{3:s}\n".format(a,b,score,info))
           else:
             out_network.write("{0} {1:f} {2}\t\t{3:s}\n".format(a,score,b,info))
      except Exception as e:
        printverbose(log,options.show,("BIANA id code error %s\n"%(e)))

       
    fd.close()


    if fileExist(options.input_nodes):
     fd=open(options.input_nodes,"r")
     for line in fd:
      if line[0] == '#':
        continue
      if options.iformat == 'netscore' :
       word=line.split()
       p=word[0]
       score=float(word[3])
       if score <= 0.0: score=options.score
       info=" ".join([str(word[i]) for i in range(4,len(word))])
      else:
       word=line.split()
       p=word[0]
       score=float(word[1])
       if score <= 0.0: score=options.score
       info=" ".join([str(word[i]) for i in range(2,len(word))])
      if new.has_key(p):
        if options.oformat ==  'sif' :
         out_nodes.write("{0} = {1}\n".format(p,";".join([x for x in new[p]])))
        else:
         for a in new[p]:
          if options.oformat == 'netscore' :
           out_nodes.write("{0}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\t\t{4:s}\n".format(a,1.,1.,score,info))
          else:
           # About STRING FORMATTING: 
           # At the left of ":"  --> 0 is the first format, 1 is the second and 2 is the third
           # At the right of ":" --> "10" is to align to the left (I have deleted it)
           #                     --> "f" is to add a float, and it indicates the number of decimal positions after the dot
           # Good examples at    --> https://pyformat.info/
           out_nodes.write("{0} {1:0.5f} {2:s}\n".format(a,score,info))
     fd.close()
    else:
     for n in nodes:
      out_nodes.write("{0}\n".format(n))

if  __name__ == "__main__":
    main()



