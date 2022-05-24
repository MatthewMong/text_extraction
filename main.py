from flask import Flask
from flask import request
from Bio import Entrez
from spacy import displacy
from spacy.matcher import Matcher
import spacy

nlp = spacy.load("en_core_web_sm")
matcher = Matcher(nlp.vocab)
matched_sents = []
pattern = [{"LIKE_NUM": True}, {"POS": "NOUN"}]
def collect_sents(matcher, doc, i, matches):
    match_id, start, end = matches[i]
    span = doc[start:end]  # Matched span
    sent = span.sent
    match_ents = [{
        "start": span.start_char - sent.start_char,
        "end": span.end_char - sent.start_char,
        "label": "MATCH",
    }]
    matched_sents.append({"text": sent.text, "ents": match_ents})

matcher.add("Subjects", [pattern], on_match=collect_sents)  # add pattern

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
            abstract = ''.join(article['Abstract']['AbstractText'])
            abstract_dict[pmid] = abstract.encode("ascii", "ignore").decode()
        else:
            without_abstract.append(pmid)
    handle.close()
    sentence_dict = {}
    for key, abstract in abstract_dict.items():
        sentence_dict[key] = []
        text = nlp(abstract)
        print(text)
        sentences = text.sents
        for match_id, start, end in matcher(text):
            print("Matched based on token shape:", text[start:end])
        for sentence in sentences:
            sentence_dict[key].append(sentence)
        
    return f'{sentence_dict}'