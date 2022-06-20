from flask import Flask
from flask import request
from Bio import Entrez
from spacy.matcher import Matcher
import spacy

nlp = spacy.load("en_core_web_trf")

n_matcher = Matcher(nlp.vocab)
sex_matcher = Matcher(nlp.vocab)
age_matcher = Matcher(nlp.vocab)
bio_fluid_matcher = Matcher(nlp.vocab)
omic_matcher = Matcher(nlp.vocab)
control_group_matcher = Matcher(nlp.vocab)
healthy_control_group_matcher = Matcher(nlp.vocab)

control = [{"LEMMA": "control"}, {"LEMMA": "group"}]
healthy_control = [{"LEMMA": "healthy"},{"LEMMA": "control"}, {"LEMMA": "group", "OP": "?"}]
control_group_matcher.add("controlGroup", [control])
healthy_control_group_matcher.add("healthyControlGroup", [healthy_control])

metabolomics = [{"LEMMA": {"IN": ["metabolomics", "metabolic", "metabolite"]}}]
proteomics = [{"LEMMA": {"IN": ["proteomic", "protein"]}}]
genomics = [{"LEMMA": {"IN": ["genomic", "genome", "gene"]}}]
microbiomics = [{"LEMMA": {"IN": ["microbiomic", "microbiome", "bacterium", "fungus", "virus"]}}]
multiomics = [{"LEMMA": {"IN": ["multiomic"]}}]

omic_matcher.add('metabolomics', [metabolomics])
omic_matcher.add('proteomics', [proteomics])
omic_matcher.add('genomics', [genomics])
omic_matcher.add('microbiomics', [microbiomics])
omic_matcher.add('multiomics', [multiomics])


bio_fluid_category = [{"LEMMA": {"IN": ["serum", "plasma", "tissue", "urine"]}}]
red_blood_cells = [{"LOWER": "red"}, {"LOWER": "blood"}, {"LEMMA": "cell"}]
rbc = [{"LOWER": "rbc"}]
bio_fluid_matcher.add('category', [bio_fluid_category])
bio_fluid_matcher.add('red_blood_cells', [red_blood_cells])
bio_fluid_matcher.add('rbc', [rbc])


age_category = [{"LEMMA":{"IN": ["adult", "teen", "child", "infant"]}}]
age_year = [{"LEMMA":{"IN": ["age, year"]}}]
age_matcher.add('category', [age_category])
age_matcher.add('year', [age_year])

subject_pattern = [{"ENT_TYPE": {"IN" : ["CARDINAL"]}, "OP": "+"}, {"POS": "ADJ", "OP": "*"}, {"POS": "NOUN", "TEXT": {"NOT_IN": ["%"]}, "LENGTH": {">": 2}}]
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
            abstract = ' '.join(article['Abstract']['AbstractText']).replace(',', '')
            abstract_dict[pmid] = abstract.encode("ascii", "ignore").decode()
        else:
            without_abstract.append(pmid)
    handle.close()
    for key, abstract in abstract_dict.items():
        text = nlp(abstract)
        potential_n = []
        matched = n_matcher(text)
        for match_id, start, end in matched:
            if len(potential_n) > 0 and end == potential_n[-1].end:
                potential_n.pop()
            potential_n.append(text[start:end])

        sexes = sex_matcher(text)
        potential_sexes = []
        for match_id, start, end in sexes:
            potential_sexes.append(text[start:end])
        if len(potential_sexes) == 0:
            potential_sexes.append("both")
        
        fluids = bio_fluid_matcher(text)
        potential_fluids = []
        for match_id, start, end in fluids:
            potential_fluids.append(text[start:end])
        
        ages = age_matcher(text)
        potential_ages = []
        for match_id, start, end in ages:
            potential_ages.append(text[start:end])
        if len(potential_ages) == 0:
            potential_ages.append("adult")
        
        omics = omic_matcher(text)
        potential_omics = []
        for match_id, start, end in omics:
            potential_omics.append(text[start:end])

        control_groups = control_group_matcher(text)
        potential_control_groups = []
        for match_id, start, end in control_groups:
            potential_control_groups.append(text[start:end])
        
        healthy_control_groups = healthy_control_group_matcher(text)
        potential_healthy_control_groups = []
        for match_id, start, end in healthy_control_groups:
            potential_healthy_control_groups.append(text[start:end])
    if (len(abstract_dict.items()) > 0):
        return f'<b>input</b><br/>{text}<br/><b>possible sample sizes:</b><br/>{potential_n}<br/><b>possible sexes:</b><br/>{set(potential_sexes)}<br/><b>possible ages:</b><br/>{set(potential_ages)}<br/><b>possible fluids:</b><br/>{set(potential_fluids)}<br/><b>possible omics:</b><br/>{set(potential_omics)}<br/><b>possible control group:</b><br/>{len(potential_control_groups) != 0}<br/><b>possible healthy control group:</b><br/>{len(potential_healthy_control_groups) != 0}'
    else:
        return '<p>no abstract</p>'


if __name__ == '__main__':
      app.run(port=8080)