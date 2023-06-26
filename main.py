import two_d_vortex_panel

if __name__ == "__main__":
    NACA_object = two_d_vortex_panel.vortex_panels("airfoils_abbreviated.json")
    NACA_object.program(0)
    NACA_object.program(1)
    NACA_object.program(2)
