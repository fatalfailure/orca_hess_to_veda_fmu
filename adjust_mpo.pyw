# -*- coding: utf-8 -*-
"""GUI launcher for Windows (.pyw).

- Runs without opening a console window.
- Shows a message box and writes a log file if an exception occurs.
"""

import os
import sys
import traceback

def _run():
    # Ensure imports work when double-clicked from Explorer
    here = os.path.abspath(os.path.dirname(__file__))
    if here and here not in sys.path:
        sys.path.insert(0, here)

    from adjust_mpo import main
    main()

if __name__ == "__main__":
    try:
        _run()
    except Exception:
        # Write a log file next to the script
        here = os.path.abspath(os.path.dirname(__file__))
        log_path = os.path.join(here, "adjust_mpo_error.log")
        tb = traceback.format_exc()
        try:
            with open(log_path, "w", encoding="utf-8") as f:
                f.write(tb)
        except Exception:
            pass

        # Show message box
        try:
            import tkinter as tk
            from tkinter import messagebox
            root = tk.Tk()
            root.withdraw()
            messagebox.showerror(
                "adjust_mpo crashed",
                "An unexpected error occurred.\n\n"
                f"A log file was written to:\n{log_path}\n\n"
                "Traceback:\n" + tb
            )
            root.destroy()
        except Exception:
            pass
