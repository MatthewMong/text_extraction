from spacy.matcher import Matcher

def get_sexes(nlp, text):
    sex_matcher = Matcher(nlp.vocab)
    female_pattern = [{"LEMMA":{"IN": ["female", "woman"]}}]
    male_pattern = [{"LEMMA":{"IN": ["male", "man"]}}]
    sex_matcher.add("female", [female_pattern])
    sex_matcher.add("male", [male_pattern])
    sexes = sex_matcher(text)
    potential_sexes = []
    for match_id, start, end in sexes:
        potential_sexes.append(text[start:end])
    if len(potential_sexes) == 0:
        potential_sexes.append("both")
    return potential_sexes