from flask import Flask
from flask import request
from Bio import Entrez
from spacy.matcher import Matcher
import spacy
from age_matcher import get_age
from control_group_matcher import get_control_groups, get_healthy_control_groups
from fluid_matcher import get_fluids
from n_matcher import get_n
from omic_matcher import get_omics
from sex_matcher import get_sexes

nlp = spacy.load("en_core_web_trf")

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
            abstract = ' '.join(article['Abstract']['AbstractText']).replace(',', '')
            abstract_dict[pmid] = abstract.encode("ascii", "ignore").decode()
        else:
            without_abstract.append(pmid)
    handle.close()
    for key, abstract in abstract_dict.items():
        text = nlp(abstract)
        potential_n = get_n(nlp, text)
        potential_sexes = get_sexes(nlp, text)
        potential_fluids = get_fluids(nlp, text)
        potential_omics = get_omics(nlp, text)
        potential_ages = get_age(nlp, text)
        potential_control_groups = get_control_groups(nlp, text)
        potential_healthy_control_groups = get_healthy_control_groups(nlp, text)
        
    if (len(abstract_dict.items()) > 0):
        return f'<b>input</b><br/>{text}<br/><b>possible sample sizes:</b><br/>{potential_n}<br/><b>possible sexes:</b><br/>{set(potential_sexes)}<br/><b>possible ages:</b><br/>{set(potential_ages)}<br/><b>possible fluids:</b><br/>{set(potential_fluids)}<br/><b>possible omics:</b><br/>{set(potential_omics)}<br/><b>possible control group:</b><br/>{len(potential_control_groups) != 0}<br/><b>possible healthy control group:</b><br/>{len(potential_healthy_control_groups) != 0}'
    else:
        return '<p>no abstract</p>'


if __name__ == '__main__':
      app.run(port=8080)