from Layer import Layer
from Multilayer import Multilayer
import pygad as pg
from numpy import arange
from Methods import *
import time
from scipy import optimize
from multiprocessing import Queue
from numba import jit



def HGGA_1(structure: Multilayer, model_psi_delta: dict, incidence_angle, queue: Queue, num_generations = 20, crossover_fraction = 0.1, sol_per_pop = None, parent_selection_type = 'rank',
                  keep_elitism = 2, mutation_type = 'random', mutation_probability = 0.5, 
                  random_mutation_min_val = -0.5, random_mutation_max_val = 0.5):

        layers_nk = structure.layers_nk
        layers = structure.layers


        missing_thicknesses = 0
        missing_thicknesses_indexes = []
        gene_space = []
        num_genes = 0
        ev_sol_per_pop = 0
        for idx, layer in enumerate(structure.layers):
            if not layer.thickness:
                missing_thicknesses += 1
                missing_thicknesses_indexes.append(idx)
                if layer.nk_file.empty:
                    if layer.roughness:
                        gene_space.append(arange(0.1, 50))
                        num_genes += 1
                        ev_sol_per_pop += ((idx+1)**2) *50
                    elif not layer.absorbing:
                        gene_space.extend((arange(1,5,0.05), arange(0,0.1,0.005), arange(0,0.01,0.0005), arange(10,1000))) #A B C thickness
                        num_genes += 4
                        ev_sol_per_pop += ((idx+1)**2) *300
                    else:
                        gene_space.extend((arange(0.1, 100), arange(0.1, 20), arange(0.1, 10, 0.5), arange(0.1, 10, 0.5), arange(0,10), arange(10,1000))) #A E0 C Eg einf thickness
                        num_genes += 6
                        ev_sol_per_pop += ((idx+1)**2)*500
                else:
                    if not layer.absorbing:
                        gene_space.append(arange(10,1000)) #A B C thickness
                        num_genes += 1
                        ev_sol_per_pop += 50
                    else:
                        gene_space.append(arange(10,1000)) #A E0 C Eg einf thickness
                        num_genes += 1
                        ev_sol_per_pop += 50
    

        if not sol_per_pop:
            sol_per_pop = ev_sol_per_pop
        else:
            sol_per_pop = sol_per_pop

        def fitness_func(ga_solution, solution, solution_idx):
            #GA_instance.mutation_probability = round(1 - GA_instance.generations_completed/num_generations, 2)
            psis_difference = 0
            deltas_difference = 0
            
            
            for idx, wvl in enumerate(model_psi_delta[case_insensitive_pick(model_psi_delta, ['wvl', 'wavelength', 'wavelength (nm)', 'nm'])]):
                list_of_layers: list[Layer] = []
                list_add = list_of_layers.append
                filled_thicknesses = 0
                used_solutions = 0
                for layer_idx, layer in enumerate(layers):
                    absorbing = layer.absorbing
                    if layer_idx in missing_thicknesses_indexes:
                        if layer.roughness:
                            layer_params = solution[used_solutions]
                            used_solutions += 1
                            if layer_idx == 0:
                                N0 = 1
                            else:
                                N0 = list_of_layers[-1].complex_refractive_indexes[1]

                            if (layer_idx+1) in missing_thicknesses_indexes:
                                if layers[layer_idx+1].absorbing:
                                    next_layer_params = solution[used_solutions:used_solutions+6]
                                    e1 = e1_TL(next_layer_params[0], next_layer_params[1], next_layer_params[2], next_layer_params[3], wavelength_to_eV(wvl), next_layer_params[4])
                                    e2 = e2_TL(next_layer_params[0], next_layer_params[1], next_layer_params[2], next_layer_params[3], wavelength_to_eV(wvl))
                                    N2 = N(n_TL(e1, e2), k_TL(e1, e2))
                                else:
                                    next_layer_params = solution[used_solutions:used_solutions+4]
                                    N2 = Cauchy_plot(next_layer_params[0], next_layer_params[1], next_layer_params[2], wavelength=wvl)

                            else:
                                N2 = N(structure.layers_nk[filled_thicknesses]['n'][idx], structure.layers_nk[filled_thicknesses]['k'][idx])
                            
                            e = EMA(eps(N=N0), eps(N=N2))
                            n = n_TL(real(e), imag(e))
                            k = k_TL(real(e), imag(e))
                            list_of_layers.append(Layer([n, k], incidence_angle, wvl, layer.absorbing, thickness=layer_params, roughness=True))
                            filled_thicknesses += 1
                        elif absorbing:
                            if layer.nk_file.empty:
                                layer_params = solution[used_solutions:used_solutions+6]
                                used_solutions += 6
                                e1 = e1_TL(layer_params[0], layer_params[1], layer_params[2], layer_params[3], wavelength_to_eV(wvl), layer_params[4])
                                e2 = e2_TL(layer_params[0], layer_params[1], layer_params[2], layer_params[3], wavelength_to_eV(wvl))
                                n = n_TL(e1, e2)
                                k = k_TL(e1, e2)
                                if n != 0 and k >= 0 and k < 20:
                                    list_add(Layer([n, k], incidence_angle, wvl, absorbing, layer.nk_file, layer_params[5]))
                                    filled_thicknesses += 1
                                else:
                                    return -1000
                            else:
                                layer_params = solution[used_solutions]
                                used_solutions += 1
                                list_add(Layer([layers_nk[filled_thicknesses]['n'][idx], layers_nk[filled_thicknesses]['k'][idx]], incidence_angle, wvl, absorbing, layer.nk_file, layer_params))
                                filled_thicknesses += 1
                        else:
                            if layer.nk_file.empty:
                                layer_params = solution[used_solutions:used_solutions+4]
                                used_solutions += 4
                                n = Cauchy_plot(layer_params[0], layer_params[1], layer_params[2], wavelength=wvl)
                                if n != 0:
                                    list_add(Layer([n, 0], incidence_angle, wvl, absorbing, layer.nk_file, layer_params[3]))
                                    filled_thicknesses+=1
                                else:
                                    return -100
                            else:
                                layer_params = solution[used_solutions]
                                used_solutions += 1
                                list_add(Layer([layers_nk[filled_thicknesses]['n'][idx], 0], incidence_angle, wvl, absorbing, layer.nk_file, layer_params))
                                filled_thicknesses += 1
                    else:
                        list_add(Layer([layers_nk[filled_thicknesses]['n'][idx], layers_nk[filled_thicknesses]['k'][idx]], incidence_angle, wvl, absorbing, layer.nk_file, layer.thickness))
                else:
                    experimental_structure = Multilayer(list_of_layers)
                    psi = experimental_structure.psi
                    delta = experimental_structure.delta
                    psis_difference += ((psi(degrees=True) - model_psi_delta[case_insensitive_pick(model_psi_delta, ['psi'])][idx]))**2
                    deltas_difference += ((delta(degrees=True) - model_psi_delta[case_insensitive_pick(model_psi_delta, ['delta'])][idx]))**2
            return 1 - sqrt((psis_difference + deltas_difference)/2)


        GA_instance = pg.GA(
            num_generations=num_generations,
            num_parents_mating=int(crossover_fraction*sol_per_pop),
            fitness_func=fitness_func,
            sol_per_pop=sol_per_pop,
            num_genes=num_genes,
            gene_space=gene_space,
            parent_selection_type=parent_selection_type,
            keep_elitism=keep_elitism,
            crossover_type='two_points',
            mutation_type=mutation_type,
            mutation_probability=mutation_probability,
            random_mutation_min_val=random_mutation_min_val,
            random_mutation_max_val=random_mutation_max_val,
            allow_duplicate_genes=False,
            stop_criteria=f'saturate_10'
        )


        start_time = time.time()
        GA_instance.run()


        solution, solution_fitness, solution_idx = GA_instance.best_solution()
        print(f'Parameters of the best solution: {solution}')
        print(f'Fitness value of the best solution: {solution_fitness}')
        print(f'It took: {round(time.time() - start_time, 2)} seconds')
        queue.put(solution)

        
def HGGA_2(structure: Multilayer, model_psi_delta: dict, incidence_angle, queue: Queue, solution = [None]):
        
        layers_nk = structure.layers_nk
        layers = structure.layers
        
        used_solutions = 0
        missing_thicknesses = 0
        missing_thicknesses_indexes = []
        bounds = []
        for idx, layer in enumerate(structure.layers):
            if not layer.thickness:
                missing_thicknesses += 1
                missing_thicknesses_indexes.append(idx)
                if layer.roughness:
                    layer_params = solution[used_solutions]
                    used_solutions += 1
                    bounds.append((0.1, layer_params+5))
                elif layer.absorbing:
                    if layer.nk_file.empty:
                        layer_params = solution[used_solutions:used_solutions+6]
                        used_solutions += 6
                        bounds.extend(((0.1,layer_params[0]+5), (0.1, layer_params[1]+2), (0.1, layer_params[2]+2), (0.1, layer_params[3]+1), (0, layer_params[4]+1), (5, layer_params[5]+100)))
                    else:
                        layer_params = solution[used_solutions]
                        used_solutions += 1
                        bounds.append((5, layer_params+100))
                else:
                    if layer.nk_file.empty:
                        layer_params = solution[used_solutions:used_solutions+4]
                        used_solutions += 4
                        bounds.extend(((0.1,layer_params[0]+1), (0, layer_params[1]+0.5), (0, layer_params[2]+0.05), (5, layer_params[3]+100)))
                    else:
                        layer_params = solution[used_solutions]
                        used_solutions += 1
                        bounds.append((5, layer_params+100))


        def fitness_func(solution):
            psis_difference = 0
            deltas_difference = 0
            
            for idx, wvl in enumerate(model_psi_delta[case_insensitive_pick(model_psi_delta, ['wvl', 'wavelength', 'wavelength (nm)', 'nm'])]):
                list_of_layers: list[Layer] = []
                filled_thicknesses = 0
                used_solutions = 0
                for layer_idx, layer in enumerate(structure.layers):
                    if layer_idx in missing_thicknesses_indexes:
                        if layer.roughness:
                            layer_params = solution[used_solutions]
                            used_solutions += 1
                            if layer_idx == 0:
                                N0 = 1
                            else:
                                N0 = list_of_layers[-1].complex_refractive_indexes[1]

                            if (layer_idx+1) in missing_thicknesses_indexes:
                                if structure.layers[layer_idx+1].absorbing:
                                    next_layer_params = solution[used_solutions:used_solutions+6]
                                    e1 = e1_TL(next_layer_params[0], next_layer_params[1], next_layer_params[2], next_layer_params[3], wavelength_to_eV(wvl), next_layer_params[4])
                                    e2 = e2_TL(next_layer_params[0], next_layer_params[1], next_layer_params[2], next_layer_params[3], wavelength_to_eV(wvl))
                                    N2 = N(n_TL(e1, e2), k_TL(e1, e2))
                                else:
                                    next_layer_params = solution[used_solutions:used_solutions+4]
                                    N2 = Cauchy_plot(next_layer_params[0], next_layer_params[1], next_layer_params[2], wavelength=wvl)
                            else:
                                N2 = N(structure.layers_nk[filled_thicknesses]['n'][idx], structure.layers_nk[filled_thicknesses]['k'][idx])
                            
                            e = EMA(eps(N=N0), eps(N=N2))
                            n = n_TL(real(e), imag(e))
                            k = k_TL(real(e), imag(e))
                            list_of_layers.append(Layer([n, k], incidence_angle, wvl, layer.absorbing, thickness=layer_params, roughness=True))
                            filled_thicknesses += 1
                        elif layer.absorbing:
                            if layer.nk_file.empty:
                                layer_params = solution[used_solutions:used_solutions+6]
                                used_solutions += 6
                                e1 = e1_TL(layer_params[0], layer_params[1], layer_params[2], layer_params[3], wavelength_to_eV(wvl), layer_params[4])
                                e2 = e2_TL(layer_params[0], layer_params[1], layer_params[2], layer_params[3], wavelength_to_eV(wvl))
                                n = n_TL(e1, e2)
                                k = k_TL(e1, e2)
                                if n != 0 and k >= 0 and k < 20:
                                    list_of_layers.append(Layer([n, k], incidence_angle, wvl, layer.absorbing, layer.nk_file, layer_params[5]))
                                    filled_thicknesses += 1
                                else:
                                    return 1000
                            else:
                                layer_params = solution[used_solutions]
                                used_solutions += 1
                                list_of_layers.append(Layer([layers_nk[filled_thicknesses]['n'][idx], layers_nk[filled_thicknesses]['k'][idx]], incidence_angle, wvl, layer.absorbing, layer.nk_file, layer_params))
                                filled_thicknesses += 1
                        else:
                            if layer.nk_file.empty:
                                layer_params = solution[used_solutions:used_solutions+4]
                                used_solutions += 4
                                n = Cauchy_plot(layer_params[0], layer_params[1], layer_params[2], wavelength=wvl)
                                if n != 0:
                                    list_of_layers.append(Layer([n, 0], incidence_angle, wvl, layer.absorbing, layer.nk_file, layer_params[3]))
                                    filled_thicknesses+=1
                                else:
                                    return 100
                            else:
                                layer_params = solution[used_solutions]
                                used_solutions += 1
                                list_of_layers.append(Layer([layers_nk[filled_thicknesses]['n'][idx], layers_nk[filled_thicknesses]['k'][idx]], incidence_angle, wvl, layer.absorbing, layer.nk_file, layer_params))
                                filled_thicknesses += 1
                    else:
                        list_of_layers.append(Layer([structure.layers_nk[filled_thicknesses]['n'][idx], structure.layers_nk[filled_thicknesses]['k'][idx]], incidence_angle, wvl, layer.absorbing, layer.nk_file, layer.thickness))
                else:
                    experimental_structure = Multilayer(list_of_layers)
                    psis_difference += ((experimental_structure.psi(degrees=True) - model_psi_delta[case_insensitive_pick(model_psi_delta, ['psi'])][idx]))**2
                    deltas_difference += ((experimental_structure.delta(degrees=True) - model_psi_delta[case_insensitive_pick(model_psi_delta, ['delta'])][idx]))**2 
            return sqrt((psis_difference + deltas_difference)/2)
        
        result = optimize.minimize(fun=fitness_func,
                                x0=solution,
                                method='L-BFGS-B',
                                bounds=bounds)
        
        
        psi_delta_dict, layers_params = create_psi_delta_dict(structure, model_psi_delta, incidence_angle, result.x)
        queue.put((layers_params, result.fun))
        queue.put(psi_delta_dict)


def create_psi_delta_dict(structure: Multilayer, model_psi_delta: dict, incidence_angle, solution):
    layers_nk = structure.layers_nk
    psi_delta = {}
    psis = []
    deltas = []
    missing_thicknesses_indexes = []
    for idx, layer in enumerate(structure.layers):
        if not layer.thickness:
            missing_thicknesses_indexes.append(idx)

    for idx, wvl in enumerate(model_psi_delta[case_insensitive_pick(model_psi_delta, ['wvl', 'wavelength', 'wavelength (nm)', 'nm'])]):
        list_of_layers: list[Layer] = []
        filled_thicknesses = 0
        used_solutions = 0
        layers_params = {}
        for layer_idx, layer in enumerate(structure.layers):
            if layer_idx in missing_thicknesses_indexes:
                if layer.roughness:
                    layer_params = solution[used_solutions]
                    used_solutions += 1
                    if layer_idx == 0:
                        N0 = 1
                    else:
                        N0 = list_of_layers[-1].complex_refractive_indexes[1]

                    if (layer_idx+1) in missing_thicknesses_indexes:
                        if structure.layers[layer_idx+1].absorbing:
                            next_layer_params = solution[used_solutions:used_solutions+6]
                            e1 = e1_TL(next_layer_params[0], next_layer_params[1], next_layer_params[2], next_layer_params[3], wavelength_to_eV(wvl), next_layer_params[4])
                            e2 = e2_TL(next_layer_params[0], next_layer_params[1], next_layer_params[2], next_layer_params[3], wavelength_to_eV(wvl))
                            N2 = N(n_TL(e1, e2), k_TL(e1, e2))
                        else:
                            next_layer_params = solution[used_solutions:used_solutions+4]
                            N2 = Cauchy_plot(next_layer_params[0], next_layer_params[1], next_layer_params[2], wavelength=wvl)
                    else:
                        N2 = N(structure.layers_nk[filled_thicknesses]['n'][idx], structure.layers_nk[filled_thicknesses]['k'][idx])
                    
                    e = EMA(eps(N=N0), eps(N=N2))
                    n = n_TL(real(e), imag(e))
                    k = k_TL(real(e), imag(e))
                    list_of_layers.append(Layer([n, k], incidence_angle, wvl, layer.absorbing, thickness=layer_params, roughness=True))
                    filled_thicknesses += 1
                elif layer.absorbing:
                    if layer.nk_file.empty:
                        layer_params = solution[used_solutions:used_solutions+6]
                        used_solutions += 6
                        e1 = e1_TL(layer_params[0], layer_params[1], layer_params[2], layer_params[3], wavelength_to_eV(wvl), layer_params[4])
                        e2 = e2_TL(layer_params[0], layer_params[1], layer_params[2], layer_params[3], wavelength_to_eV(wvl))
                        n = n_TL(e1, e2)
                        k = k_TL(e1, e2)
                        list_of_layers.append(Layer([n, k], incidence_angle, wvl, layer.absorbing, layer.nk_file, layer_params[5]))
                        filled_thicknesses+=1
                    else:
                        layer_params = solution[used_solutions]
                        used_solutions += 1
                        list_of_layers.append(Layer([layers_nk[filled_thicknesses]['n'][idx], layers_nk[filled_thicknesses]['k'][idx]], incidence_angle, wvl, layer.absorbing, layer.nk_file, layer_params))
                        filled_thicknesses += 1


                else:
                    if layer.nk_file.empty:
                        layer_params = solution[used_solutions:used_solutions+4]
                        used_solutions += 4
                        n = Cauchy_plot(layer_params[0], layer_params[1], layer_params[2], wavelength=wvl)
                        list_of_layers.append(Layer([n, 0], incidence_angle, wvl, layer.absorbing, layer.nk_file, layer_params[3]))
                        filled_thicknesses+=1
                    else:
                        layer_params = solution[used_solutions]
                        used_solutions += 1
                        list_of_layers.append(Layer([layers_nk[filled_thicknesses]['n'][idx], 0], incidence_angle, wvl, layer.absorbing, layer.nk_file, layer_params))
                        filled_thicknesses += 1



                layers_params[layer_idx] = layer_params

            else:
                list_of_layers.append(Layer([structure.layers_nk[filled_thicknesses]['n'][idx], structure.layers_nk[filled_thicknesses]['k'][idx]], incidence_angle, wvl, layer.absorbing, layer.nk_file, layer.thickness))
        else:
            experimental_structure = Multilayer(list_of_layers)
            psis.append(experimental_structure.psi(degrees=True))
            deltas.append(experimental_structure.delta(degrees=True))
    else:
        psi_delta['wvl'] = model_psi_delta[case_insensitive_pick(model_psi_delta, ['wvl', 'wavelength', 'wavelength (nm)', 'nm'])]
        psi_delta['psi'] = psis
        psi_delta['delta'] = deltas
        return psi_delta, layers_params