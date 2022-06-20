from spacy.matcher import Matcher

def get_fluids(nlp, text):
    bio_fluid_matcher = Matcher(nlp.vocab)
    bio_fluid_category = [{"LEMMA": {"IN": ["serum", "plasma", "tissue", "urine"]}}]
    red_blood_cells = [{"LOWER": "red"}, {"LOWER": "blood"}, {"LEMMA": "cell"}]
    rbc = [{"LOWER": "rbc"}]
    bio_fluid_matcher.add('category', [bio_fluid_category])
    bio_fluid_matcher.add('red_blood_cells', [red_blood_cells])
    bio_fluid_matcher.add('rbc', [rbc])
    fluids = bio_fluid_matcher(text)
    potential_fluids = []
    for match_id, start, end in fluids:
        potential_fluids.append(text[start:end])
    return potential_fluids