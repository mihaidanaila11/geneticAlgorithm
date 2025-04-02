import tkinter as tk

class windowInterface:
    def __init__(self, genAlgorithm = None):
        def startHandler(event):
            arguments = []
            for input in self.inputs:
                inputValue = input.get()

                if not inputValue:
                    print("Invalid arguments!")
                    return
                try:
                    arguments.append(int(inputValue))
                except:
                    print("Invalid arguments!")
                    return
            genAlgorithm(*arguments).evolve()

        def createInputWidget(labelText):
            widget = tk.Label(text=labelText)
            entry = tk.Entry()
            return [widget, entry]

        self.window = tk.Tk()

        widgets = []

        widgets += createInputWidget("Population Size")
        widgets += createInputWidget("Left limit of the interval")
        widgets += createInputWidget("Right limit of the interval")
        widgets += [tk.Label(text="Function parameters f(x)=a*x^2 + b*x + c")]
        widgets += createInputWidget("a")
        widgets += createInputWidget("b")
        widgets += createInputWidget("c")
        widgets += createInputWidget("Precision")
        widgets += createInputWidget("Crossover Probability")
        widgets += createInputWidget("Mutation Probability")
        widgets += createInputWidget("Number of generations")

        self.inputs = [widget for widget in widgets if isinstance(widget, tk.Entry)]


        startButton = tk.Button(text="Start Algorithm")
        startButton.bind("<Button-1>", startHandler)
        widgets.append(startButton)

        for widget in widgets:
            widget.pack()

        self.widgets = widgets

        self.window.mainloop() 