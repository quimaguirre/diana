import sys, os

"""
    NetworkAnalysis
    2017 Joaquim Aguirre-Plans 
    Structural Bioinformatics Laboratory
    Universitat Pompeu Fabra
"""

def translate(input_network, input_nodes, translation_file, input_format, output_format, output_network, output_nodes):

    if not fileExist(translation_file):
        print("Translation file is missing\n")
        sys.exit(10)

    new={}

    # Fill new with the translations
    if fileExist(translation_file):
        ft=open(translation_file,"r")
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


    if fileExist(input_network):
        fd=open(input_network,"r")
        out_network=output_network
        if not output_network == sys.stdout: out_network=open(output_network,"w")

        # seeds=set()
        # if fileExist(options.seed):
        #  fn=open(options.seed,'r')
        #  for line in fn:
        #   word=line.split()
        #   seeds.add(word[0])
        #  fn.close()

        nodes=set()

        for line in fd:

            # Process the input network
            skip=False
            if input_format == 'multi-fields' :
                word=line.strip().split("\t")
                p=word[0] # Node 1
                q=word[1] # Node 2
                try:
                    score=word[2]       # Code for Janet # This is not the score, it is the source, but it is ok
                except:
                    score=1.0
            elif input_format == 'sif':
                word=line.split()
                if len(word)>2:
                    p=word[0]
                    q=word[2]
                    score=word[1]
                else:skip=True
            elif input_format == 'raw':
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
            info="\t".join([str(word[i]) for i in range(3,len(word))])       # Code for Janet
            if not translation:
                new.setdefault(p,set()).add(p)
                new.setdefault(q,set()).add(q)

            # Print the output network
            for a in new[p]:
                nodes.add(a)
                for b in new[q]:
                    nodes.add(b)
                    if output_format == 'multi-fields' :
                        out_network.write("{0}\t{1}\t{2}\t{3}\n".format(a,b,score,info))       # Code for Janet
                    elif output_format == 'sif' :
                        if not skip: out_network.write("{0}\t{1}\t{2}\n".format(a,score,b))
                    else:
                        out_network.write("{0} {1:f} {2}\t\t{3:s}\n".format(a,score,b,info))
        fd.close()
        if not output_network == sys.stdout: out_network.close()




    # Process the input nodes
    if fileExist(input_nodes):

        out_network=output_nodes
        if not output_nodes == sys.stdout: out_network=open(output_nodes,"w")
        if output_format != 'multi-fields':
            out_network.write("#Id Score P-value\n")

        fd=open(input_nodes,"r")

        for line in fd:
            if line[0] == '#':
                continue
            if input_format == 'multi-fields' :
                word=line.split()
                p=word[0]
                info="\t".join([str(word[i]) for i in range(1,len(word))])
            else:
                word=line.split()
                p=word[0]
                score=float(word[1])
                info=" ".join([str(word[i]) for i in range(2,len(word))])

            # Print the output nodes
            if new.has_key(p):
                if output_format ==  'sif' :
                    out_network.write("{0} = {1}\n".format(p,";".join([x for x in new[p]])))
                else:
                    for a in new[p]:
                        if output_format == 'multi-fields':
                            out_network.write("{0}\t{1}\n".format(a,info))
                        else:
                            # About STRING FORMATTING: 
                            # At the left of ":"  --> 0 is the first format, 1 is the second and 2 is the third
                            # At the right of ":" --> "10" is to align to the left (I have deleted it)
                            #                     --> "f" is to add a float, and it indicates the number of decimal positions after the dot
                            # Good examples at    --> https://pyformat.info/
                            out_network.write("{0} {1:f} {2:s}\n".format(a,score,info))
        fd.close()

    # else:
    #     for a in nodes:
    #         if a in seeds: score=1.0
    #         else:          score=options.score
    #         if output_format == 'netscore' :
    #             if a is not None: out_network.write("{0}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\n".format(a,1.,1.,score))
    #         elif output_format ==  'sif' :
    #             if translation:
    #                 if new.has_key(a): 
    #                     if a is not None: out_network.write("{0} = {1}\n".format(a,";".join([x for x in new[a]])))
    #                 else:
    #                     if a is not None: out_network.write("{0}\n".format(a))
    #             else:
    #                 if a is not None: out_network.write("{0}\n".format(a))
    #         else:
    #             if a is not None: out_network.write("{0} {1:10.5f}\n".format(a,score))

    out_network.close()

    return


def fileExist (file):               #Checks if a file exists AND is a file
    if file is not None: return os.path.exists(file) and os.path.isfile(file)
    else: return False
