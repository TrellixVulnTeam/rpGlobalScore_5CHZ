import tarfile
import rpSBML
import libsbml

def rpExtractJSONfromSBML(inputTar, pathway_id='rp_pathway'):
    json_out = {}
    with tarfile.open(inputTar, mode='r:xz') as in_tf:
        for member in in_tf.getmembers():
            if not member.name=='':
                member_name = member.name.replace('.xml', '').replace('.rpsbml', '').replace('.sbml', '')
                rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(in_tf.extractfile(member).read().decode('utf-8')))
                json_out[member_name] = {}
                #Groups
                groups = rpsbml.model.getPlugin('groups')
                rp_pathway = groups.getGroup(pathway_id)
                json_out[member_name][pathway_id] = rpsbml.readBRSYNTHAnnotation(rp_pathway.getAnnotation())
                json_out[member_name]['reactions'] = {}
                for reac_mem in rp_pathway.getListOfMembers():
                    reac = rpsbml.model.getReaction(reac_mem.getIdRef())
                    json_out[member_name]['reactions'][reac_mem.getIdRef()] = rpsbml.readBRSYNTHAnnotation(reac.getAnnotation())
                #FBC
                fbc = rpsbml.model.getPlugin('fbc')
                json_out[member_name]['objectives'] = {}
                for objective in fbc.getListOfObjectives():
                   json_out[member_name]['objectives'][objective.getId()] = rpsbml.readBRSYNTHAnnotation(objective.getAnnotation())
    return json_out
