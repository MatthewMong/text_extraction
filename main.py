from flask import Flask
from flask import request
from Bio import Entrez
from spacy.matcher import Matcher
from spacy import displacy
import spacy

nlp = spacy.load("en_core_web_sm")

n_matcher = Matcher(nlp.vocab)
sex_matcher = Matcher(nlp.vocab)
# subject_pattern = [{"LIKE_NUM": True}, {"POS": "ADJ", "OP": "*"}, {"POS": "NOUN", "DEP": {"IN": ["acomp","pobj", "ROOT", "nsubjpass", "nsubj", "conj","compound"]}}]

subject_pattern = [{"LIKE_NUM": True}, {"POS": "ADJ", "OP": "*"}, {"POS": "NOUN", "TEXT": {"NOT_IN": ["%"]}}]
acronym_pattern = [{"LIKE_NUM": True}, {"IS_UPPER": True}]
n_pattern = [{"TEXT": {"REGEX": "^n="}}]
n_spaces_pattern = [{"LOWER": "n"}, {"TEXT": "="}, {"LIKE_NUM": True}]
n_matcher.add("subjects", [subject_pattern])  # add pattern
n_matcher.add("n_string", [n_pattern])  # add pattern
n_matcher.add("n", [n_spaces_pattern])  # add pattern
n_matcher.add("acronym", [acronym_pattern])  # add pattern

female_pattern = [{"LEMMA":{"IN": ["female", "woman"]}}]
male_pattern = [{"LEMMA":{"IN": ["male", "man"]}}]
sex_matcher.add("female", [female_pattern])
sex_matcher.add("male", [male_pattern])

email = ''
app = Flask(__name__)

@app.route('/paper')
def papers():
    doi = request.args.get('doi')
    Entrez.email = email
    handle = Entrez.elink(dbfrom="pubmed", db="pmc", linkname="pubmed_pmc", id=doi, retmode="xml")
    id_return = Entrez.read(handle)
    handle.close()
    handle = Entrez.efetch(db="pmc", id=id_return[0]['LinkSetDb'][0]['Link'][0]['Id']) 
    records = handle.read()
    handle.close()
    return f'{records}'

@app.route("/")
def entrance():
    doi = request.args.getlist('doi')
    abstract_dict = {}
    without_abstract = []
    Entrez.email = email
    handle = Entrez.efetch(db="pubmed", id=''.join(doi),rettype="xml", retmode="text") 
    records = Entrez.read(handle)
    for pubmed_article in records['PubmedArticle']:
        pmid = int(str(pubmed_article['MedlineCitation']['PMID']))
        article = pubmed_article['MedlineCitation']['Article']
        if 'Abstract' in article:
            abstract = ' '.join(article['Abstract']['AbstractText'])
            abstract_dict[pmid] = abstract.encode("ascii", "ignore").decode()
        else:
            without_abstract.append(pmid)
    handle.close()
    for key, abstract in abstract_dict.items():
        text = nlp(abstract)
        # for tok in text:
        #     print(tok.i, tok, tok.pos_, tok.dep_, tok.head.i, sep="\t")
        potential_n = []
        matched = n_matcher(text)
        for match_id, start, end in matched:
            potential_n.append(text[start:end])
        sexes = sex_matcher(text)
        potential_sexes = []
        for match_id, start, end in sexes:
            potential_sexes.append(text[start:end])
    return f'{text}</br>possible sample sizes: </br>{potential_n}</br>possible sexes: </br>{set(potential_sexes)}'


if __name__ == '__main__':
      app.run(port=8080)