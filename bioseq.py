#!/usr/bin/python3
#######################################################################
# bioseq.py: Módulo para leitura e gravação de arquivos fasta         #
# Autor: Felipe Nunes de Freitas Lima                                 #
# Data: 02/10/2017                                                    #
#######################################################################

def leFasta(file):
	try:
		seq = {}
		f = open(file, 'r')
		lines = f.readlines()
		f.close()

		for line in lines:
			if line[0] == '>':
				lineCurrent = lines.index(line) + 1
				seq[line.strip()] = ''
				while lineCurrent < len(lines) and lines[lineCurrent][0] != '>':
					seq[line.strip()] += lines[lineCurrent].strip()
					lineCurrent += 1
		return seq
	except FileNotFoundError:
		print ("O arquivo especificado não existe!")

def gravaFasta(file, seq):
	if type(seq) == dict:
		f = open(file,'w')
		contentFile = ''

		for key in seq.keys():

			if key[0] != '>':
				contentFile += '>'
			contentFile += key + '\n'
				
			countColumns = 1
			for column in seq[key]:
				contentFile += column
				countColumns += 1				

				if countColumns == 60:
					countColumns = 1
					contentFile += '\n'
			contentFile += '\n\n'

		f.write(contentFile)
		f.close
	else:
		raise ValueError('O parâmetro seq precisa ser um dicionário. Foi passado como parâmetro um tipo '+str(type(seq))+'.')

def lePHD(file):
	try:
		filePHD = open(file,'r')
		seq = list()
		seqFasta = dict()
		seqName = ''
		beginDNA = False

		for line in filePHD.readlines():
			if line[0:14] == 'BEGIN_SEQUENCE':
				seqName = line[15:-1]
			elif line[0:9] == 'BEGIN_DNA':
				beginDNA = True
			elif line[0:7] == 'END_DNA':
				break
			elif beginDNA == True:
				seqToAppend = line.split(' ')
				seqToAppend.pop()
				seq.append(seqToAppend)
		
		filePHD.close()

		fileQuality = open(seqName+'.qual','w')
		fileQuality.write('>'+seqName+'\n')
		
		seqFasta[seqName] = ''

		for line in seq:
			seqFasta[seqName] += line[0]
			fileQuality.write(line[1]+'\n')

		fileQuality.close()
		gravaFasta(seqName+'.fasta', seqFasta)
			
	except FileNotFoundError:
		print("O arquivo específicado não existe!")
