import customtkinter as tk
from tkinter import filedialog as fd
from typing import List
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from Thickness_get_methods import *
from multiprocessing import Queue, freeze_support
from threading import Thread
from queue import Empty

queue = Queue()


class FirstWindow(tk.CTkFrame):
    def __init__(self, parent):
        super().__init__(parent)

        self.result = None
        self.file_to_transition = None
        self.experimentalStructure = None
        self.transitionButton = None
        self.layers = None
        self.canvas = None
        self.ax = None
        self.fig = None
        self.psi_delta_plot_frame = None
        self.selected_psi_delta = None
        self.filename = None

        psiDeltaFrame = tk.CTkFrame(self)

        psiDeltaLabel = tk.CTkLabel(psiDeltaFrame,
                                    text='Select file containing Psi and Delta values')
        psiDeltaLabel.pack()

        self.psiDeltaButton = tk.CTkButton(psiDeltaFrame,
                                           text='Select file',
                                           command=self.select_psiDelta_file)
        self.psiDeltaButton.pack()

        psiDeltaFrame.grid(column=2, row=0, sticky=tk.NE, padx=(200, 5), pady=(5, 200))

        # Entry for incidence angle
        incidenceAngleFrame = tk.CTkFrame(self)

        angleLabel = tk.CTkLabel(incidenceAngleFrame,
                                 text='Incidence angle(degrees)')
        angleLabel.grid(column=0, row=0, padx=5, pady=2)

        self.angleEntry = tk.CTkEntry(incidenceAngleFrame,
                                      textvariable=tk.StringVar(value=70),
                                      width=28)
        self.angleEntry.grid(column=0, row=1, padx=5, pady=2)

        wavelengthLabel = tk.CTkLabel(incidenceAngleFrame,
                                      text='Wavelength(nm)')
        wavelengthLabel.grid(column=1, row=0, padx=5, pady=2)

        self.wavelengthEntry = tk.CTkEntry(incidenceAngleFrame,
                                           textvariable=tk.StringVar(value=589.3),
                                           width=45)
        self.wavelengthEntry.grid(column=1, row=1, padx=5, pady=2)

        incidenceAngleFrame.grid(column=0, row=0, sticky=tk.NW, padx=5, pady=5)

        # Set of buttons
        # Add layer
        layerButtonsFrame = tk.CTkFrame(self)

        self.addLayerButton = tk.CTkButton(layerButtonsFrame,
                                           text='Add layer',
                                           command=self.addLayer)
        self.addLayerButton.pack(side='left')

        self.addRoughnessButton = tk.CTkButton(layerButtonsFrame,
                                               text='Add roughness',
                                               command=self.addRoughness)
        self.addRoughnessButton.pack(side='left')

        layerButtonsFrame.grid(column=0, row=1, padx=5, pady=5, sticky=tk.SW)

        # Create layers frame with substrate layer
        self.layersFrame = tk.CTkScrollableFrame(self, width=290, height=250)

        self.layersFrames: List[LayerFrame] = []

        substrateLayer = LayerFrame(self.layersFrame,
                                    'Substrate')
        substrateLayer.grid(row=0, column=0)

        self.layersFrames.append(substrateLayer)
        self.layersFrame.grid(column=0, row=2, sticky=tk.NW, padx=(2, 20), pady=(2, 10))

        self.optionsMenuFrame = tk.CTkFrame(self)

        self.optionsMenu = tk.CTkOptionMenu(self.optionsMenuFrame,
                                            values=['Find thickness', 'Generate Psi and Delta',
                                                    'Delta range transition'],
                                            command=self.optionsMenuCallback)
        self.optionsMenu.grid(column=0, row=0)
        self.optionsMenu.set('Find thickness')

        self.runButton = tk.CTkButton(self.optionsMenuFrame, text='Run', command=self.runButtonCommand)
        self.runButton.grid(column=0, row=1)

        self.optionsMenuFrame.grid(column=2, row=1, sticky=tk.E, padx=5, pady=5)

        self.resultsFrame = tk.CTkFrame(self)

        self.resultsLabelList: list[list[tk.CTkLabel]] = []


        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=2)
        self.grid_columnconfigure(2, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=3)

        self.pack(fill='both', expand=True)

    # Adds layer on top
    def addLayer(self):
        layer = LayerFrame(self.layersFrame)
        self.layersFrames.insert(0, layer)

        for layer_idx, layer in enumerate(self.layersFrames):
            layer.grid(row=layer_idx)

    def addRoughness(self):
        roughnessLayer = LayerFrame(self.layersFrame,
                                    name='Roughness')
        roughnessLayer.nFrame.pack_forget()
        roughnessLayer.kFrame.pack_forget()
        roughnessLayer.absorbFrame.pack_forget()
        roughnessLayer.nk_buttonFrame.pack_forget()
        self.layersFrames.insert(0, roughnessLayer)

        for layer_idx, layer in enumerate(self.layersFrames):
            layer.grid(row=layer_idx)

    def select_psiDelta_file(self):
        filetypes = [('Comma-Separated Values', '*.csv')]

        self.filename = fd.askopenfilename(
            title='Open a csv file',
            initialdir='/',
            filetypes=filetypes)

        if self.psiDeltaButton.cget('text') != 'Select file':
            self.create_base_psi_delta_plot(self.canvas, self.ax)
        else:
            self.create_psi_delta_plot()

            self.create_base_psi_delta_plot(self.canvas, self.ax)

        self.selected_psi_delta = pd.read_csv(self.filename, sep=';')

        self.psiDeltaButton.configure(text=self.filename.split('/')[len(self.filename.split('/')) - 1])

    def create_psi_delta_plot(self):

        self.psi_delta_plot_frame = tk.CTkFrame(self)
        self.fig = Figure(figsize=(5, 5),
                          dpi=100)

        self.ax = self.fig.add_axes([0.1, 0.1, 0.8, 0.8])

        self.canvas = FigureCanvasTkAgg(self.fig,
                                        master=self.psi_delta_plot_frame)
        self.canvas.get_tk_widget().pack()
        self.canvas.draw()

        toolbar = NavigationToolbar2Tk(self.canvas,
                                       self.psi_delta_plot_frame)
        toolbar.update()
        self.psi_delta_plot_frame.grid(column=1, row=2, padx=10, pady=10)

    def create_base_psi_delta_plot(self, canvas, ax):

        file = pd.read_csv(self.filename, sep=';')

        wavelengthValues = file[case_insensitive_pick(file, ['wvl', 'wavelength', 'wavelength (nm)', 'nm'])]
        psiValues = file[case_insensitive_pick(file, ['psi'])]
        deltaValues = file[case_insensitive_pick(file, ['delta'])]

        ax.clear()
        ax.plot(wavelengthValues, psiValues, color='red', label='Psi')
        ax.plot(wavelengthValues, deltaValues, color='green', label='Delta')
        ax.legend()
        ax.grid()
        canvas.draw()

    def update_psi_delta_plot(self, canvas, ax, psi_delta: dict):
        self.create_base_psi_delta_plot(canvas, ax)
        for idx, wvl in enumerate(psi_delta[case_insensitive_pick(psi_delta, ['wvl', 'wavelength', 'wavelength (nm)', 'nm'])]):
            ax.scatter(wvl, psi_delta[case_insensitive_pick(psi_delta, ['psi'])][idx], c='tab:red', marker='s')
            ax.scatter(wvl, psi_delta[case_insensitive_pick(psi_delta, ['delta'])][idx], c='tab:green', marker='^')
        canvas.draw()

    def get_values(self):
        self.layers: List[Layer] = []

        for layerFrameIdx, layerFrame in enumerate(self.layersFrames):
            nk_file = layerFrame.selected_file
            roughness = False
            try:
                if layerFrame.nEntry.get() != '':
                    n = float(layerFrame.nEntry.get())
                else:
                    n = 1

                if layerFrame.kEntry.get() != '':
                    k = float(layerFrame.kEntry.get())
                else:
                    k = 1

                if layerFrame.absorbCheckVar.get() == 0:
                    absorbing = False
                else:
                    absorbing = True
            except:
                continue

            if layerFrame.name == 'Substrate':
                thickness = 1
                absorbing = True
            elif layerFrame.name == 'Roughness':
                roughness = True
                absorbing = False
                if layerFrame.thicknessEntry.get():
                    thickness = float(layerFrame.thicknessEntry.get())
                else:
                    thickness = None

                try:
                    previousLayer = self.layersFrames[layerFrameIdx-1]
                    N0 = N(float(previousLayer.nEntry.get()), float(previousLayer.kEntry.get()))
                except:
                    N0 = 1
                    
                nextLayer = self.layersFrames[layerFrameIdx+1]
                try:
                    N2 = N(float(nextLayer.nEntry.get()), float(nextLayer.kEntry.get()))
                except:
                    N2 = 1
                e = EMA(N0**2, N2**2)
                n = n_TL(real(e), imag(e))
                k = k_TL(real(e), imag(e))
            elif layerFrame.thicknessEntry.get() != '':
                thickness = float(layerFrame.thicknessEntry.get())
            else:
                self.layers.append(Layer([n, 0],
                                        float(self.angleEntry.get()),
                                        float(self.wavelengthEntry.get()),
                                        absorbing,
                                        nk_file))
                continue


            self.layers.append(Layer([n, k],
                                     float(self.angleEntry.get()),
                                     float(self.wavelengthEntry.get()),
                                     absorbing,
                                     nk_file,
                                     thickness,
                                     roughness))

        else:
            self.experimentalStructure = Multilayer(self.layers)

    def optionsMenuCallback(self, choice):
        optionsMenuValue = self.optionsMenu.get()
        if optionsMenuValue == 'Delta range transition':
            self.runButton.grid(row=2)
            self.transitionButton = tk.CTkButton(self.optionsMenuFrame, text='Select file',
                                                 command=self.select_file_to_transition)
            self.transitionButton.grid(column=0, row=1)
            self.optionsMenuFrame.update()
        else:
            try:
                self.runButton.grid(row=1)
                self.transitionButton.destroy()
            except:
                pass

    def runButtonCommand(self):
        optionsMenuValue = self.optionsMenu.get()
        if optionsMenuValue == 'Find thickness':
            self.findThickness()
        elif optionsMenuValue == 'Generate Psi and Delta':
            self.generate_psi_delta()
        elif optionsMenuValue == 'Delta range transition':
            self.file_transition()

    def findThickness(self):
        while not queue.empty():
            try:
                queue.get(block=False)
            except Empty:
                continue
        self.get_values()

        HGGA_process_1 = Thread(name='HGGA_1', target=HGGA_1, args=(self.experimentalStructure,
                                                                    self.selected_psi_delta,
                                                                    float(self.angleEntry.get()),
                                                                    queue,
                                                                    30,))
        HGGA_process_1.start()
        HGGA_process_1.join()
        result1 = queue.get()

        HGGA_process_2 = Thread(name='HGGA_2', target=HGGA_2, args=(self.experimentalStructure,
                                                                    self.selected_psi_delta,
                                                                    float(self.angleEntry.get()),
                                                                    queue,
                                                                    result1,))
        HGGA_process_2.start()
        HGGA_process_2.join()
        self.layers_params, self.result = queue.get()
        self.psi_delta = queue.get()
        print(f'layers_params: {self.layers_params}')
        self.update_psi_delta_plot(self.canvas, self.ax, self.psi_delta)
        self.update_results(thickness=True)

    def generate_psi_delta(self):
        self.get_values()

        self.update_results(psi_delta=True)

    def select_file_to_transition(self):
        filetypes = [('Comma-Separated Values', '*.csv')]

        self.file_to_transition = fd.askopenfilename(
            title='Open a csv file',
            initialdir='/',
            filetypes=filetypes)

        self.transitionButton.configure(
            text=self.file_to_transition.split('/')[len(self.file_to_transition.split('/')) - 1])

    def file_transition(self):
        transitioned_file = CompleteEASE_to_normal(self.file_to_transition, sep=';')
        transitioned_file.to_csv(f'{self.file_to_transition.split(".csv")[0]}_transitioned.csv', sep=';')

    def update_results(self, thickness: bool = False, psi_delta: bool = False):
        self.resultsFrame.grid(column=2, row=2)
        try:
            for widget in self.resultsFrame.winfo_children():
                widget.destroy()
        except:
            pass
        if thickness:
            idx = 0
            for layer_idx in self.layers_params.keys():
                tk.CTkLabel(self.resultsFrame, text=f'Parameters of layer number {layer_idx}:').grid(column=0, row=idx, padx=3, pady=2)
                try:
                    tk.CTkLabel(self.resultsFrame, text=f'{"thickness"}: {round(self.layers_params[layer_idx], 4)}').grid(column=0, row=idx+1, padx=3, pady=2)
                except:
                    if len(self.layers_params[layer_idx]) == 4:
                        tk.CTkLabel(self.resultsFrame, text=f'{"A"}: {round(self.layers_params[layer_idx][0], 4)}').grid(column=0, row=idx+1, padx=3, pady=2)
                        tk.CTkLabel(self.resultsFrame, text=f'{"B"}: {round(self.layers_params[layer_idx][1], 4)}').grid(column=1, row=idx+1, padx=3, pady=2)
                        tk.CTkLabel(self.resultsFrame, text=f'{"C"}: {round(self.layers_params[layer_idx][2], 4)}').grid(column=2, row=idx+1, padx=3, pady=2)
                        tk.CTkLabel(self.resultsFrame, text=f'{"thickness"}: {round(self.layers_params[layer_idx][3], 4)}').grid(column=3, row=idx+1, padx=3, pady=2)
                    else:
                        tk.CTkLabel(self.resultsFrame, text=f'{"A"}: {round(self.layers_params[layer_idx][0], 4)}').grid(column=0, row=idx+1, padx=3, pady=2)
                        tk.CTkLabel(self.resultsFrame, text=f'{"E0"}: {round(self.layers_params[layer_idx][1], 4)}').grid(column=1, row=idx+1, padx=3, pady=2)
                        tk.CTkLabel(self.resultsFrame, text=f'{"C"}: {round(self.layers_params[layer_idx][2], 4)}').grid(column=2, row=idx+1, padx=3, pady=2)
                        tk.CTkLabel(self.resultsFrame, text=f'{"Eg"}: {round(self.layers_params[layer_idx][3], 4)}').grid(column=3, row=idx+1, padx=3, pady=2)
                        tk.CTkLabel(self.resultsFrame, text=f'{"e_inf"}: {round(self.layers_params[layer_idx][4], 4)}').grid(column=4, row=idx+1, padx=3, pady=2)
                        tk.CTkLabel(self.resultsFrame, text=f'{"thickness"}: {round(self.layers_params[layer_idx][5], 4)}').grid(column=5, row=idx+1, padx=3, pady=2)
                idx += 2
            else:
                cohesion = round(((1-self.result)/1)*100, 2)
                if cohesion < 0:
                    cohesion = 0
                tk.CTkLabel(self.resultsFrame, text=f'Result: {round(self.result, 6)}(cohesion: {cohesion}%)').grid(column=0, row=idx, padx=3, pady=2)
                self.resultsFrame.update()
        elif psi_delta:
            tk.CTkLabel(self.resultsFrame, text=f'Psi: {self.experimentalStructure.psi(degrees=True)}').grid(column=0, row=0, padx=3, pady=2)
            tk.CTkLabel(self.resultsFrame, text=f'Delta: {self.experimentalStructure.delta(degrees=True)}').grid(column=0, row=1, padx=3, pady=2)
            self.resultsFrame.update()
        else:
            raise ValueError




class LayerFrame(tk.CTkFrame):
    def __init__(self, parent, name=None):
        # Set name
        self.name = name
        self.selected_file = pd.DataFrame(None)
        super().__init__(parent)
        if name:
            tk.CTkLabel(self, text=name).pack()
        else:
            name = tk.CTkEntry(self, textvariable=tk.StringVar(), width=6)
            name.pack()
            self.name = name.get()


        # Button for complex refractive index selection
        self.nk_buttonFrame = tk.CTkFrame(self)
        self.nk_button = tk.CTkButton(self.nk_buttonFrame,
                                      text='Select N data',
                                      command=self.select_nk_file)
        self.nk_button.pack()
        self.nk_buttonFrame.pack()

        # Thickness entry
        thicknessFrame = tk.CTkFrame(self)
        thicknessLabel = tk.CTkLabel(thicknessFrame,
                                     text='Thickness(nm)').pack()
        self.thicknessEntry = tk.CTkEntry(thicknessFrame,
                                          textvariable=tk.StringVar(),
                                          width=50)
        self.thicknessEntry.pack()
        thicknessFrame.pack(side='left')

        # n entry
        self.nFrame = tk.CTkFrame(self)
        nLabel = tk.CTkLabel(self.nFrame,
                             text='n').pack()
        self.nEntry = tk.CTkEntry(self.nFrame,
                                  textvariable=tk.StringVar(),
                                  width=50)
        self.nEntry.pack()
        self.nFrame.pack(side='left')

        # k entry
        self.kFrame = tk.CTkFrame(self)
        kLabel = tk.CTkLabel(self.kFrame,
                             text='k').pack()
        self.kEntry = tk.CTkEntry(self.kFrame,
                                  textvariable=tk.StringVar(),
                                  width=50)
        self.kEntry.pack()
        self.kFrame.pack(side='left')

        # absorbing button
        self.absorbFrame = tk.CTkFrame(self)
        absorbLabel = tk.CTkLabel(self.absorbFrame,
                                  text="Absorbing").grid(column=0, row=0)
        self.absorbCheckVar = tk.IntVar()  #check box variable
        self.absorbCheckButton = tk.CTkCheckBox(self.absorbFrame,
                                                variable=self.absorbCheckVar,
                                                # width=8,
                                                onvalue=1,
                                                offvalue=0,
                                                text='')
        self.absorbCheckButton.grid(column=0, row=1)
        self.absorbFrame.pack(side='left')

    # command for nk file selection
    def select_nk_file(self):
        filetypes = [('Comma-Separated Values', '*.csv')]

        filename = fd.askopenfilename(
            title='Open a csv file',
            initialdir='/',
            filetypes=filetypes)

        self.nk_button.configure(text=filename.split('/')[len(filename.split('/')) - 1])

        self.selected_file = pd.read_csv(filename, sep=';')

   


class MainWindow:
    def __init__(self, master):
        mainframe = tk.CTkFrame(master)
        mainframe.pack(fill='both', expand=True)
        self.windowNum = 0

        self.framelist: List[tk.CTkFrame] = []
        self.framelist.append(FirstWindow(mainframe))
        for frameIdx, frame in enumerate(self.framelist):
            if frameIdx != self.windowNum:
                frame.forget()


if __name__ == '__main__':
    freeze_support()
    root = tk.CTk()
    root.geometry('1920x1000')
    root.state('zoomed')
    root.title('WSE 2023')
    window = MainWindow(root)
    root.mainloop()
