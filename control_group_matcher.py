from spacy.matcher import Matcher

def get_control_groups(nlp, text):
    control_group_matcher = Matcher(nlp.vocab)
    control = [{"LEMMA": "control"},{"POS": "ADJ", "OP": "*"}, {"LEMMA": {"IN":["group", "subject"]}}]
    controls = [{"LOWER": "controls", "POS": "NOUN"}]
    control_group_matcher.add("controlGroup", [control, controls])
    control_groups = control_group_matcher(text)
    potential_control_groups = []
    for match_id, start, end in control_groups:
        potential_control_groups.append(text[start:end])
    return potential_control_groups

def get_healthy_control_groups(nlp, text):
    healthy_control_group_matcher = Matcher(nlp.vocab)
    healthy_control = [{"LEMMA": "healthy"},{"POS": "ADJ", "OP": "*"},{"LEMMA": "control", "POS": "NOUN"}, {"LEMMA": {"IN": ["group", "subject"]}, "OP": "?"}]
    healthy_control_group_matcher.add("healthyControlGroup", [healthy_control])
    healthy_control_groups = healthy_control_group_matcher(text)
    potential_healthy_control_groups = []
    for match_id, start, end in healthy_control_groups:
        potential_healthy_control_groups.append(text[start:end])
    return potential_healthy_control_groups
