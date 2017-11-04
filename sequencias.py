#!/usr/bin/python3
#######################################################################
# sequencias.py: Classes para leitura de sequências biológicas	      #
# Autor: Felipe Nunes de Freitas Lima                                 #
# Data: 27/10/2017                                                    #
#######################################################################
import bioseq

class Sequencia:

	def __init__(self,cadeia = '', genBank = '', descricao = ''):
		self.cadeia = cadeia
		self.genBank = genBank
		self.descricao = descricao

	def salvarFasta(self,caminho):
		key = "> "
		if not self.cadeia:
			raise Exception("Não é possível salvar o arquivo pois a cadeia está vazia")
		if self.genBank:
			key += "GenBank: "+self.genBank
		if self.descricao:
			key += " Descrição: "+self.descricao
		seq = {key : self.cadeia}
		bioseq.gravaFasta(caminho, seq)

	def abrirFasta(self,caminho):
		try:
			seq = bioseq.leFasta(caminho)
		except:
			raise Exception("Erro abrir arquivo!")

		if len(seq) == 1:

			cabecalho = list(seq.keys())[0]

			self.cadeia = list(seq.values())[0]
			if '|' in cabecalho:
				self.genBank = cabecalho.split('|')[0]
				for descricao in cabecalho.split('|'):
					if descricao != self.genBank:
						self.descricao += descricao
			else:
				self.descricao = list(seq.keys())[0]
		else:
			raise Exception("O arquivo FASTA deve possuir apenas uma sequência!")

	def getTamanho(self):
		return len(self.cadeia)

	def getComposicaoAbsoluta(self):
		aminoacids = {}
		for aminoacid in self.cadeia:
			if aminoacid in aminoacids.keys():
				aminoacids[aminoacid] += 1
			else:
				aminoacids[aminoacid] = 1

		return aminoacids

	def getComposicaoRelativa(self):
		aminoacids = {}
		aminoacidsRelative = {}

		for aminoacid in self.cadeia:
			if aminoacid in aminoacids.keys():
				aminoacids[aminoacid] += 1
			else:
				aminoacids[aminoacid] = 1

		for aminoacid in aminoacids:
			aminoacidsRelative[aminoacid] = aminoacids[aminoacid] / len(aminoacids)

		return aminoacidsRelative

class RNA(Sequencia):

	def __init__(self,cadeia = '', genBank = '', descricao = ''):
		Sequencia.__init__(self, cadeia, genBank, descricao)

		self.dicTradutor = {  'TGT':'C', 'TGC':'C', 'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
						 'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V', 'GTC':'V',
						 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCA':'S', 'TCG':'S', 'TCC':'S', 'CCT':'P', 'CCC':'P',
						 'CCG':'P', 'CCA':'P', 'ACT':'T', 'ACG':'T', 'ACA':'T', 'ACC':'T', 'GCT':'A', 'GCC':'A',
						 'GCG':'A', 'GCA':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*', 'CAT':'H', 'CAC':'H',
						 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D',
						 'GAA':'E', 'GAG':'E', 'TGA':'*', 'TGG':'W', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
						 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

	def getProteina(self):
		proteinas = list()

		for frame in range(0,3):
			i = frame
			proteinas.append(Proteina())
			while i + 2 < len(self.cadeia):
				try:
					proteinas[len(proteinas) - 1].cadeia += self.dicTradutor[self.cadeia[i:i+3]]
				except KeyError:
					print('Uma trinca de nucleotídeos não possui correspondência de aminoácido. Verifique se os dados da cadeia estão corretos')
					return
				i += 3

		for frame in range(0,3):
			i = frame
			proteinas.append(Proteina())
			while i + 2 < len(self.cadeia):
				try:
					proteinas[len(proteinas) - 1].cadeia += self.dicTradutor[self.cadeia[::-1][i:i+3]]
				except KeyError:
					print('\n\nUma trinca de nucleotídeos não possui correspondência de aminoácido. Verifique se os dados da cadeia estão corretos')
					return
				i += 3

		return proteinas

class Proteina(Sequencia):

	def __init__(self,cadeia = '', genBank = '', descricao = ''):
		Sequencia.__init__(self, cadeia, genBank, descricao)

class DNA(Sequencia):

	def __init__(self,cadeia = '', genBank = '', descricao = ''):
		Sequencia.__init__(self, cadeia, genBank, descricao)

if __name__ == "__main__":
	print("Testes:\n")
	rna = RNA()
	rna.abrirFasta('seq.fas')
	print("Abrindo seq.fas: \n")
	print("Lendo todas as proteínas possíveis geradas a partir dessa sequência:")
	for proteina in rna.getProteina():
		print("\n"+proteina.cadeia)

	print("\nSalvando proteína e abrindo o arquivo em outra instância: \n")
	proteina = rna.getProteina()[0]
	proteina.salvarFasta('proteinaTeste.fas')
	proteinaTeste = Proteina()
	proteinaTeste.abrirFasta('proteinaTeste.fas')
	print(proteinaTeste.cadeia)
	print("Tamanho: "+str(proteinaTeste.getTamanho()))
	print("Frequência Absoluta: "+str(proteinaTeste.getComposicaoAbsoluta()))
	print("Frequência Relativa: "+str(proteinaTeste.getComposicaoRelativa()))
