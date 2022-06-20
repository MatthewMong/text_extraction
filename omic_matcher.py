from spacy.matcher import Matcher

def get_omics(nlp, text):
    omic_matcher = Matcher(nlp.vocab)
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
 
    omics = omic_matcher(text)
    potential_omics = []
    for match_id, start, end in omics:
        potential_omics.append(text[start:end])
    return potential_omics