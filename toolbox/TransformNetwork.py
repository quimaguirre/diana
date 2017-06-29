import sys
import argparse
import os


def main():

    options = parse_user_arguments()
    translate(options)


def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Translate a bPPI network into your selected codes",
        epilog      = "@oliva's lab 2014")
    parser.add_argument('-i','--network_input_file',dest='input_edge',action = 'store',default='edges.txt',
                        help = 'Input file of network (default is edges.txt)')
    parser.add_argument('-n','--nodes_input_file',dest='input_node',action = 'store',default='nodes.txt',
                        help = 'Input file of nodes (default is nodes.txt)')
    parser.add_argument('-trans','--translation_of_nodes_file',dest='translation_file',action = 'store',default=None,
                        help = 'File with the translation of codes from BIANA to the selected type for all nodes')
    parser.add_argument('-iformat','--input_format',dest='iformat',action = 'store',default='guild',
                        help = 'Format of input files of edges:\tguild (default) \tnetscore \tsif\traw')
    parser.add_argument('-oformat','--output_format',dest='oformat',action = 'store',default='guild',
                        help = 'Format of output files of edges:\tguild (default) \tnetscore \tsif')
    parser.add_argument('-oe','--output_file_edges',dest='output_edge',action = 'store', default=sys.stdout,
                        help = 'Output file with edges(default is standard output)')
    parser.add_argument('-on','--output_file_nodes',dest='output_node',action = 'store', default=sys.stdout,
                        help = 'Output file with nodes(default is standard output)')
    parser.add_argument('-score','--Non-seed_score',dest='score',action = 'store',default='0.0',type=float,
                        help = 'Score of non-seed nodes (default is 0.0)')
    parser.add_argument('-seed','--seed_file',dest='seed',action = 'store',default=None,
                        help = 'File with seeds. It will output in the file of nodes with scores 1.0')
    options=parser.parse_args()

    return options

def fileExist (file):               #Checks if a file exists AND is a file
 if file is not None: return os.path.exists(file) and os.path.isfile(file)
 else: return False

def translate(options):
    if not fileExist(options.input_edge):
     print("File with input network is missing\n")
     sys.exit(10)

    new={}
    if fileExist(options.translation_file):
     ft=open(options.translation_file,"r")
     for line in ft:
      word=line.strip().split("\t")
      node=word[0]
      new.setdefault(node,set())
      for transet in word[1].split("'"):
       for trans in transet.split(","):
        if trans != "" and trans != ",":
          name="_".join([str(x) for x in trans.split()])
          new[node].add(name)
     ft.close()
     translation=True
    else:
     translation=False

    fd=open(options.input_edge,"r")
    out_network=options.output_edge
    if not options.output_edge == sys.stdout: out_network=open(options.output_edge,"w")

    seeds=set()
    if fileExist(options.seed):
     fn=open(options.seed,'r')
     for line in fn:
      word=line.split()
      seeds.add(word[0])
     fn.close()

    nodes=set()
    for line in fd:
      skip=False
      if options.iformat == 'netscore' :
       word=line.split()
       p=word[0]
       q=word[1]
       try:
        score=float(word[2])
        #score=word[2]       # Code for Janet
       except:
        score=1.0
      elif options.iformat == 'sif':
       word=line.split()
       if len(word)>2:
        p=word[0]
        q=word[2]
        score=1.0
       else:skip=True
      elif options.iformat == 'raw':
       word=line.split()
       if len(word)>1:
        p=word[0]
        q=word[1]
        score=1.0
       else:skip=True
      else:
       word=line.split()
       p=word[0]
       score=float(word[1])
       q=word[2]
      info="".join([str(word[i]) for i in range(3,len(word))])
      #info=" ".join([str(word[i]) for i in range(3,len(word))])       # Code for Janet
      if not translation:
       new.setdefault(p,set()).add(p)
       new.setdefault(q,set()).add(q)
      for a in new[p]:
       nodes.add(a)
       for b in new[q]:
          nodes.add(b)
          if options.oformat == 'netscore' :
            out_network.write("\t{0}\t{1}\t{2:f}\t\t{3:s}\n".format(a,b,score,info))
            #out_network.write("{0}\t{1}\t{2}\t{3}\n".format(a,b,score,info))       # Code for Janet
          elif options.oformat == 'sif' :
            if not skip: out_network.write("{0} interaction {1}\n".format(a,b))
          else:
            out_network.write("{0} {1:f} {2}\t\t{3:s}\n".format(a,score,b,info))
    fd.close()
    if not options.output_edge == sys.stdout: out_network.close()



    out_network=options.output_node
    if not options.output_node == sys.stdout: out_network=open(options.output_node,"w")
    if options.oformat ==  'sif' : out_network.write("#%s\n"%(options.translation_file))


    if fileExist(options.input_node):
     fd=open(options.input_node,"r")
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
         out_network.write("{0} = {1}\n".format(p,";".join([x for x in new[p]])))
        else:
         for a in new[p]:
          if options.oformat == 'netscore' :
           out_network.write("{0}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\t\t{4:s}\n".format(a,1.,1.,score,info))
          else:
           # About STRING FORMATTING: 
           # At the left of ":"  --> 0 is the first format, 1 is the second and 2 is the third
           # At the right of ":" --> "10" is to align to the left (I have deleted it)
           #                     --> "f" is to add a float, and it indicates the number of decimal positions after the dot
           # Good examples at    --> https://pyformat.info/
           out_network.write("{0} {1:0.5f} {2:s}\n".format(a,score,info))
     fd.close()
    else:
     for a in nodes:
          if a in seeds: score=1.0
          else:          score=options.score
          if options.oformat == 'netscore' :
           if a is not None: out_network.write("{0}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\n".format(a,1.,1.,score))
          elif options.oformat ==  'sif' :
            if translation:
                if new.has_key(a): 
                  if a is not None: out_network.write("{0} = {1}\n".format(a,";".join([x for x in new[a]])))
                else:
                  if a is not None: out_network.write("{0}\n".format(a))
            else:
                if a is not None: out_network.write("{0}\n".format(a))
          else:
           if a is not None: out_network.write("{0} {1:10.5f}\n".format(a,score))

    out_network.close()

if  __name__ == "__main__":
    main()



