from spacy.matcher import Matcher

def get_age(nlp, text):
    age_matcher = Matcher(nlp.vocab)
    age_category = [{"LEMMA":{"IN": ["adult", "teen", "child", "infant"]}}]
    age_year = [{"LEMMA":{"IN": ["age, year"]}}]
    age_matcher.add('category', [age_category])
    age_matcher.add('year', [age_year])
    ages = age_matcher(text)
    potential_ages = []
    for match_id, start, end in ages:
        potential_ages.append(text[start:end])
    if len(potential_ages) == 0:
        potential_ages.append("adult")
    return potential_ages