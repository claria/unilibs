#!/usr/bin/env python

import sys
import os
import numpy
import glob
import argparse
import StringIO

def main():
    """
    Create LHgrid file from HERAFitter lhapdf.block.txt
    """
    parser = argparse.ArgumentParser(description='Create LHgrid File.')
    parser.add_argument("-i", "--input_folder", help="Input Folder")
    parser.add_argument("-o", "--output_filepath", help="Output Folder", default='./test.LHgrid')
    parser.add_argument("-l", "--limit_members", help="Limit Members Folder", default=0, type=int)
    parser.add_argument("--eig", help="Produce eigenvector pdf", action='store_true')

    args = vars(parser.parse_args())

    if args['eig']:
        produce_eig_lhgrid(**args)
    else:
        produce_general_lhgrid(**args)

def produce_eig_lhgrid(**args):

    member_folder = args['input_folder'] + '/MEMBER_00'
    alphas = query(member_folder + '/minuit.out.txt','alphas')[-1].split()[2]
    order = query(member_folder + '/steering.txt','Order')[0].split()[2].replace('\'','')
    npar = query(member_folder + '/minuit.out.txt','NPAR')[-1].split()[6]

    lhgrids = []
    lhgrids.append(args['input_folder'] + '/MEMBER_00'+ '/lhapdf.block.txt')
    lhgrids += sorted(glob.glob(args['input_folder'] + '/MEMBER_00'+ '/pdfs_*.lhgrid'))
    print lhgrids

    header= get_header(nmembers=len(lhgrids)-1,
                       order=order,
                       alphas_values=len(lhgrids)*(alphas,))

    with open(args['output_filepath'], 'w') as f:
        f.write(header)
        for n, lhgrid in enumerate(lhgrids):
            print lhgrid
            griddef, pdfblock = get_lhapdf_block(lhgrid)
            if n == 0:
                griddef.seek(0)
                f.write(griddef.read())
            pdfblock.seek(0)
            f.write(pdfblock.read())
        f.write("'End:'\n")



def produce_general_lhgrid(**args):

    member_folders = sorted(glob.glob(args['input_folder'] + '/*'))

    if args['limit_members'] != 0:
        member_folders = member_folders[0:args['limit_members']]

    metadata = {'alphas' : [], 'order' : [], 'npar' : []}
    valid_pdfs = []

    #Get all the physics parameters
    print "Begin checking PDF folders and extracting metadata"

    for member_folder in member_folders:
        try:
            alphas = query(member_folder + '/minuit.out.txt','alphas')[-1].split()[2]
            order = query(member_folder + '/steering.txt','Order')[0].split()[2].replace('\'','')
            npar = query(member_folder + '/minuit.out.txt','NPAR')[-1].split()[6]
            #Check if lhapdf.block.txt exists
            with open(member_folder + '/lhapdf.block.txt') as file:
                pass  #Potential for unhandled exception
        except:
            print "Not a valid lhgrid folder", member_folder
            continue
        #PDF seems to be valid
        valid_pdfs.append(member_folder)
        metadata['alphas'].append(alphas)
        metadata['order'].append(order)
        metadata['npar'].append(npar)

    print "Scanned all folders"
    print "There are {0} valid pdfs".format(len(valid_pdfs))


    header= get_header(nmembers=len(valid_pdfs)-1,
                       order=metadata['order'][0],
                       alphas_values=metadata['alphas'])

    print "This is the used header file."
    print "Check carefully."

    print header


    print "Begin writing the lhgrid blocks to file"

    with open(args['output_filepath'], 'w') as f:
        f.write(header)
        for n, member_folder in enumerate(valid_pdfs):
            print "Member {0}: {1}".format(n, member_folder)
            print "Member {0}: AlphasMz: {0}, Order: {0}".format(n, metadata['alphas'][n], metadata['order'][n])
            griddef, pdfblock = get_lhapdf_block(member_folder + '/lhapdf.block.txt')
            if n == 0:
                griddef.seek(0)
                f.write(griddef.read())
            pdfblock.seek(0)
            f.write(pdfblock.read())
        f.write("'End:'\n")



def query(file, pattern):
    result = []
    with open(file, 'r') as f:
        for line in f:
            if pattern in line:
                result.append(line)
    return result

def get_header(nmembers=0, order='nlo', alphas_values=None):

    header="""'Version' '5.8'
'Description: '
'HERAPDF fit: varaint'
'member 0 - central'
' '
'Alphas:'
'Variable','{order}','EvolCode'
1,91.187,1.40,4.75,180.0
'MinMax:'
{nmembers},1
1.E-06,1.,1.0,200000000.
'QCDparams:'
{nmembers},1
0.338,0.243
'Parameterlist:'
'list',{nmembers},1
{alphas_values}
'Evolution:'
'{order}',1.9,1.0
'HERAGRID10'
{nmembers},0
"""
    header =  header.format(nmembers=nmembers,alphas_values='\n'.join(alphas_values),order=order.lower())
    return header


def make_lhgrid(input_files, output_file, alphasmz, evolution):

    q2valpdf = numpy.zeros(161)
    xvalpdf = numpy.zeros(161)
    alphas = numpy.zeros(161)
    buffer = numpy.zeros((8,161,161))

    if os.path.isfile(output_file):
        print "Warning. File exists"
        os.remove(output_file)

    nmember = len(input_files)-1

    header="""'Version' '5.8'
'Description: '
'HERAPDF fit: varaint'
'member 0 - central'
' '
'Alphas:'
'Variable','{evolution}','EvolCode'
1,91.187,1.40,4.75,180.0
'MinMax:'
{nmember},1
1.E-06,1.,1.0,200000000.
'QCDparams:'
{nmember},1
0.338,0.243
'Parameterlist:'
'list',{nmember},1
{alphasmz}
'Evolution:'
'{evolution}',1.9,1.0
'HERAGRID10'
{nmember},0
"""
    header =  header.format(nmember=nmember,alphasmz='\n'.join((nmember+1)*(alphasmz,)),evolution=evolution.lower())

    print "Using the following header. Please check carefully!"
    print header


    with open(output_file, 'a') as f:
        f.write(header)

        with open(output_file,'a') as f:
            if n == 0:
                for iq2 in range(0,23):
                    line = []
                    for j in range(0,7):
                        line.append(q2valpdf[iq2*7+j])
                    f.write(''.join('{0:15.7E}'.format(x) for x in line) + '\n')
                for jx in range(0,23):
                    line = []
                    for j in range(0,7):
                        line.append(xvalpdf[jx*7+j])
                    f.write(''.join('{0:15.7E}'.format(x) for x in line) + '\n')
                for iq2 in range(0,23):
                    line = []
                    for j in range(0,7):
                        line.append(alphas[iq2*7+j])
                    f.write(''.join('{0:15.7E}'.format(x) for x in line) + '\n')
            for k in range(0,8):
                for i in range(0,161):
                        for j in range(0,23):
                            line = []
                            for ii in range(0,7):
                                line.append(buffer[k][j*7 + ii][i])
                            f.write(''.join('{0:15.7E}'.format(x) for x in line) + '\n')

    with open(output_file,'a') as f:
        f.write("'End:'\n")

def get_lhapdf_block(filepath):

    q2valpdf = numpy.zeros(161)
    xvalpdf = numpy.zeros(161)
    alphas = numpy.zeros(161)
    buffer = numpy.zeros((8,161,161))

    with open(filepath) as f:
        for iq2 in range(0,23):
            line = f.readline().split()
            for j in range(0,7):
                q2valpdf[iq2*7 + j] = line[j]
        for jx in range(0,23):
            line = f.readline().split()
            for j in range(0,7):
                xvalpdf[(jx)*7+j] = line[j]
        for iq2 in range(0,23):
            line = f.readline().split()
            for j in range(0,7):
                alphas[(iq2)*7+j] = line[j]
        for i in range(0,161):
            for j in range(0,161):
                line = f.readline().split()
                for k in range(0,8):
                    buffer[k][j][i] = line[k]

    griddef = StringIO.StringIO()
    pdfblock = StringIO.StringIO()
    for iq2 in range(0,23):
        line = []
        for j in range(0,7):
            line.append(q2valpdf[iq2*7+j])
        griddef.write(''.join('{0:15.7E}'.format(x) for x in line) + '\n')
    for jx in range(0,23):
        line = []
        for j in range(0,7):
            line.append(xvalpdf[jx*7+j])
        griddef.write(''.join('{0:15.7E}'.format(x) for x in line) + '\n')
    for iq2 in range(0,23):
        line = []
        for j in range(0,7):
            line.append(alphas[iq2*7+j])
        griddef.write(''.join('{0:15.7E}'.format(x) for x in line) + '\n')

    for k in range(0,8):
        for i in range(0,161):
                for j in range(0,23):
                    line = []
                    for ii in range(0,7):
                        line.append(buffer[k][j*7 + ii][i])
                    pdfblock.write(''.join('{0:15.7E}'.format(x) for x in line) + '\n')

    return griddef, pdfblock

def write_blocks(blocklist, filepath='test.LHgrid'):

    with open(output_file, 'a') as f:
        f.write(header)

        with open(output_file,'a') as f:
            if n == 0:
                for iq2 in range(0,23):
                    line = []
                    for j in range(0,7):
                        line.append(q2valpdf[iq2*7+j])
                    f.write(''.join('{0:15.7E}'.format(x) for x in line) + '\n')
                for jx in range(0,23):
                    line = []
                    for j in range(0,7):
                        line.append(xvalpdf[jx*7+j])
                    f.write(''.join('{0:15.7E}'.format(x) for x in line) + '\n')
                for iq2 in range(0,23):
                    line = []
                    for j in range(0,7):
                        line.append(alphas[iq2*7+j])
                    f.write(''.join('{0:15.7E}'.format(x) for x in line) + '\n')
            for k in range(0,8):
                for i in range(0,161):
                        for j in range(0,23):
                            line = []
                            for ii in range(0,7):
                                line.append(buffer[k][j*7 + ii][i])
                            f.write(''.join('{0:15.7E}'.format(x) for x in line) + '\n')

    with open(output_file,'a') as f:
        f.write("'End:'\n")
    pass

if __name__ == '__main__':
    sys.exit(main())
