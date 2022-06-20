from spacy.matcher import Matcher

def get_control_groups(nlp, text):
    control_group_matcher = Matcher(nlp.vocab)
    control = [{"LEMMA": "control"}, {"LEMMA": "group", "OP": "?"}]
    control_group_matcher.add("controlGroup", [control])
    control_groups = control_group_matcher(text)
    potential_control_groups = []
    for match_id, start, end in control_groups:
        potential_control_groups.append(text[start:end])
    return potential_control_groups

def get_healthy_control_groups(nlp, text):
    healthy_control_group_matcher = Matcher(nlp.vocab)
    healthy_control = [{"LEMMA": "healthy"},{"LEMMA": "control"}, {"LEMMA": "group", "OP": "?"}]
    healthy_control_group_matcher.add("healthyControlGroup", [healthy_control])
    healthy_control_groups = healthy_control_group_matcher(text)
    potential_healthy_control_groups = []
    for match_id, start, end in healthy_control_groups:
        potential_healthy_control_groups.append(text[start:end])
    return potential_healthy_control_groups