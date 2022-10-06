def update_namelist(conf, nml, section="", list_update=True):
    if section:
        print("update_namelist called")
        for key in nml[section].keys():
            if hasattr(conf, "%s"%key.upper()):
                nml[section][key][0] = float(getattr(conf, key.upper()))
                print(section, key, 'updated')
            for k in range(len(nml[section][key])):
                    if hasattr(conf, "%s_%d"%(key.upper(), k)):
                        nml[section][key][k] = float(getattr(conf, "%s_%d"%(key.upper(), k)))
                        print(section, key, k, 'updated')
            
    else:
        for section in nml.keys():
            for key in nml[section].keys():
                for k in range(len(nml[section][key])):
                    if hasattr(conf, "%s_%s_%d"%(section.upper(), key.upper(), k)):
                        nml[section][key][k] = float(getattr(conf, "%s_%s_%d"%(section.upper(), key.upper(), k)))
                        print(section, key, k, 'updated')

