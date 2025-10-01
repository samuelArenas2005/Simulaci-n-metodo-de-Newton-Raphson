"""
Simple Python interpreter GUI using tkinter.

Features:
- Input area to type Python code
- Run button to execute the code in a safe-ish local namespace
- Output area showing stdout and stderr
- Clear buttons for input and output
- Basic execution in a background thread to keep the UI responsive

NOTES & SAFETY:
- This runs arbitrary Python code (exec). Do NOT run untrusted code.
- It's meant for learning / experimenting locally.

How to run:
$ python tkinter_python_interpreter.py

"""
import tkinter as tk
from tkinter import ttk, font, messagebox
import threading
import io
import sys
import traceback
from contextlib import redirect_stdout, redirect_stderr


class SimpleInterpreter(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Python intérprete - Tkinter")
        self.geometry("900x600")
        self._create_widgets()
        # Namespace for exec (persist between runs)
        self.user_globals = {"__name__": "__main__"}

    def _create_widgets(self):
        # Fonts
        monospace = font.Font(family="Consolas" if "Consolas" in font.families() else "Courier", size=11)

        # Top frame for controls
        top = ttk.Frame(self)
        top.pack(side=tk.TOP, fill=tk.X, padx=6, pady=6)

        run_btn = ttk.Button(top, text="Ejecutar (F5)", command=self._on_run)
        run_btn.pack(side=tk.LEFT)

        clear_in_btn = ttk.Button(top, text="Limpiar entrada", command=lambda: self.input_text.delete("1.0", tk.END))
        clear_in_btn.pack(side=tk.LEFT, padx=(6, 0))

        clear_out_btn = ttk.Button(top, text="Limpiar salida", command=lambda: self.output_text.delete("1.0", tk.END))
        clear_out_btn.pack(side=tk.LEFT, padx=(6, 0))

        status_lbl = ttk.Label(top, text="Intérprete simple — cuidado: ejecuta código arbitrario.")
        status_lbl.pack(side=tk.RIGHT)

        # Paned window to split input/output
        paned = ttk.PanedWindow(self, orient=tk.VERTICAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=6, pady=(0,6))

        # Input area
        input_frame = ttk.LabelFrame(paned, text="Código (entrada)")
        self.input_text = tk.Text(input_frame, wrap=tk.NONE, height=15, font=monospace, undo=True)
        self.input_text.pack(fill=tk.BOTH, expand=True)

        # Line numbers for input (simple)
        self._attach_line_numbers(input_frame)

        paned.add(input_frame, weight=3)

        # Output area
        output_frame = ttk.LabelFrame(paned, text="Salida (stdout / stderr)")
        self.output_text = tk.Text(output_frame, wrap=tk.NONE, height=12, font=monospace, state=tk.NORMAL)
        self.output_text.pack(fill=tk.BOTH, expand=True)

        paned.add(output_frame, weight=2)

        # Bindings
        self.bind_all('<F5>', lambda e: self._on_run())
        self.input_text.bind('<Control-Return>', lambda e: self._on_run())

    def _attach_line_numbers(self, parent_frame):
        # A small line number implementation on the left of the input_text
        ln_canvas = tk.Text(parent_frame, width=4, padx=4, takefocus=0, border=0, background="#f0f0f0", state=tk.DISABLED)
        ln_canvas.pack(side=tk.LEFT, fill=tk.Y)  # Cambio: usar pack en lugar de place
        
        # Agregar padding izquierdo al input_text
        self.input_text.configure(padx=8)  # Agregar padding
        
        def update_line_numbers(event=None):
            # compute number of lines and show
            ln_canvas.configure(state=tk.NORMAL)
            ln_canvas.delete('1.0', tk.END)
            lines = int(self.input_text.index('end-1c').split('.')[0])
            ln_str = '\n'.join(str(i) for i in range(1, lines + 1)) + '\n'
            ln_canvas.insert('1.0', ln_str)
            ln_canvas.configure(state=tk.DISABLED)

            # sync scrolling
            ln_canvas.yview_moveto(self.input_text.yview()[0])

        # attach scroll sync
        self.input_text.config(yscrollcommand=lambda *args: (ln_canvas.yview_moveto(args[0])) )
        self.input_text.bind('<KeyRelease>', update_line_numbers)
        self.input_text.bind('<MouseWheel>', update_line_numbers)
        self.input_text.bind('<Button-1>', update_line_numbers)
        self.input_text.bind('<ButtonRelease-1>', update_line_numbers)
        update_line_numbers()

    def _on_run(self):
        code = self.input_text.get('1.0', 'end-1c')
        if not code.strip():
            return

        # disable run button to prevent double runs
        # start background thread
        thread = threading.Thread(target=self._execute_code, args=(code,), daemon=True)
        thread.start()

    def _execute_code(self, code):
        # Capture stdout and stderr
        stdout_buf = io.StringIO()
        stderr_buf = io.StringIO()
        try:
            with redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
                # compile first to show syntax errors nicely
                compiled = compile(code, '<user_code>', 'exec')
                exec(compiled, self.user_globals)
        except Exception:
            # write traceback to stderr buffer
            traceback.print_exc(file=stderr_buf)

        out = stdout_buf.getvalue()
        err = stderr_buf.getvalue()

        # Append result to output widget on the main thread
        self.after(0, self._append_output, out, err)

    def _append_output(self, out, err):
        if out:
            self.output_text.insert(tk.END, out)
        if err:
            self.output_text.insert(tk.END, err)
        # Auto-scroll to end
        self.output_text.see(tk.END)


if __name__ == '__main__':
    app = SimpleInterpreter()
    app.mainloop()
