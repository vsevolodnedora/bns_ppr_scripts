# contains set of methods to load/dump select/sort/get files

from lists import *

class MakePath:

    def __init__(self):
        pass

    @staticmethod
    def outflow(sim, outflow = '_0'):
        return Paths.ppr_sims + sim + '/' + 'outflow' + outflow + '/'

    @staticmethod
    def collated(sim):
        return Paths.ppr_sims + sim + '/' + 'collated' + '/'

    @staticmethod
    def collated2(sim):
        return Paths.ppr_sims + sim + '/' + 'collated2' + '/'

    @staticmethod
    def waveforms(sim):
        return Paths.ppr_sims + sim + '/' + 'waveforms' + '/'

def get_profiles(sim, time_list = [], it_list = [], n_more=0,
                 ftype = '.cylh5', time_units='s', add_path='profiles/'):
    '''
    Looks at all the available profiles in the sim.dir + '/profiles/'.
    Using the .it_time_file, interpolates [times] for all [iterations]

    '''

    from general import interpoate_it_form_time, interpoate_time_form_it, find_nearest_index
    import numpy as np
    from glob import glob

    if len(time_list) > 0 and len(it_list) > 0: raise NameError("Only time_list or it_list can be specified")

    files = sorted(glob(Paths.ppr_sims + sim + "/" + add_path + "*{}".format(ftype)),
                   key=lambda x: int(x.split('/')[-1].replace("{}".format(ftype), "")))
    iterations = [int(it.split('/')[-1].replace("{}".format(ftype), "")) for it in files]
    get_times = interpoate_time_form_it(iterations, Paths.ppr_sims + sim + '/' + Files.it_time, time_units)

    print('|------------------------------------------------------------------------|')
    print("ALL TIMES: \n{}".format(get_times))
    print("ALL ITERATIONS: \n{}".format(np.array(iterations, dtype=int)))
    print('|------------------------------------------------------------------------|')

    # if nothing specified, return all the files found and all the timestpes
    if not any(time_list) and not any(it_list):
        return get_times, files

    # if times to find are given, but not iterations, - compute iterations for those times and find files
    elif any(time_list) and not any(it_list):
        get_iterations = interpoate_it_form_time(time_list, Paths.ppr_sims + sim + '/' + Files.it_time, time_units)

    # if times are not given, but iterations to use are, - compute times for those iterations
    elif not any(time_list) and any(it_list):
        if np.array(it_list, dtype=int).min() < np.array(iterations, dtype=int).min():
            raise ValueError("Given it_list it:{} < iterations.min()"
                             .format(np.array(it_list, dtype=int).min(), np.array(iterations, dtype=int).min()))
        if np.array(it_list, dtype=int).max() > np.array(iterations, dtype=int).max():
            raise ValueError("Given it_list it:{} < iterations.min()"
                             .format(np.array(it_list, dtype=int).max(), np.array(iterations, dtype=int).max()))
        get_iterations = it_list
    else:
        raise IOError("Input is not recongized.")

    # select those iterations that are present in the 'iterations' list of files.
    available_iterations = []
    files_of_available_iterations = []
    if len(get_iterations) < 1:
        raise ValueError("No iterations for times:{} were found".format(time_list))

    for it in get_iterations:
        available_iterations.append(iterations[find_nearest_index(np.array(iterations, dtype=int), it)])
        fname = Paths.ppr_sims + sim + "/" + add_path + "{}{}".format(available_iterations[-1], ftype)
        files_of_available_iterations.append(fname)

    if len(available_iterations) < 1:
        raise ValueError("No available iterations seleted from required list: {}".format(get_iterations))
    if len(available_iterations) < len(get_iterations):
        raise ValueError("N of available it:{} != N of required:{}".format(len(available_iterations), len(get_iterations)))

    # if just one time/iteration is given and n_more > 0, get MORE iterations (consecutive)
    if n_more > 0:
        if len(time_list) == 1 or len(it_list) == 1:
            idx = find_nearest_index(np.array(iterations, dtype=int), int(available_iterations[0]))
            for i in range(1, n_more):
                fname = Paths.ppr_sims + sim + "/" + add_path + "{}{}".format(iterations[int(idx + i)], ftype)
                available_iterations.append(iterations[idx + i])
                files_of_available_iterations.append(fname)

    available_times = interpoate_time_form_it(available_iterations,
                                              Paths.ppr_sims + sim + '/' + Files.it_time, time_units)

    return available_times, files_of_available_iterations