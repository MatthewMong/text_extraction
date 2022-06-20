from spacy.matcher import Matcher

def get_n(nlp, text):
    n_matcher = Matcher(nlp.vocab)
    subject_pattern = [{"ENT_TYPE": {"IN" : ["CARDINAL"]}, "OP": "+"}, {"POS": "ADJ", "OP": "*"}, {"POS": "NOUN", "TEXT": {"NOT_IN": ["%"]}, "LENGTH": {">": 2}}]
    acronym_pattern = [{"LIKE_NUM": True}, {"IS_UPPER": True}]
    n_pattern = [{"TEXT": {"REGEX": "^n="}}]
    n_spaces_pattern = [{"LOWER": "n"}, {"TEXT": "="}, {"LIKE_NUM": True}]
    n_matcher.add("subjects", [subject_pattern])  # add pattern
    n_matcher.add("n_string", [n_pattern])  # add pattern
    n_matcher.add("n", [n_spaces_pattern])  # add pattern
    n_matcher.add("acronym", [acronym_pattern])  # add pattern
    potential_n = []
    matched = n_matcher(text)
    for match_id, start, end in matched:
        if len(potential_n) > 0 and end == potential_n[-1].end:
            potential_n.pop()
        potential_n.append(text[start:end])
    return potential_n