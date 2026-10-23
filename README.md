# Automated SMLM microscope

This repository describes a custom-built SMLM extension for a Nikon Ti-E microscope, designed for high-throughput TIRF-SMLM measurements. An additional RCM 1 unit is attached to the microscope’s left port for rapid sample screening.

The setup is controlled using NIS-Elements, extended with custom [NIS-Elements macros](https://github.com/Rickmic/Receptor_Quantification/tree/main/high_throughput_microscopy/NIS_macro). Three Arduino Mega microcontrollers provide additional hardware control ([Arduino scripts](https://github.com/Rickmic/Receptor_Quantification/tree/main/high_throughput_microscopy/Arduino)):

- **Arduino 1:** Controls the TIRF translation stage and triggers UV laser pulses.
- **Arduino 2:** Controls the DPSS lasers via RS-232 and positions a servo-driven neutral-density filter wheel.
- **Arduino 3:** Controls the custom-built piezo slider system within the RCM unit.

All three Arduinos interface with NIS-Elements through a National Instruments DAQ card using digital and analog signals.

<p align="center">
 <b>SMLM beam path of the high-throughput microscope</b>
</p>

![SMLM beam path of the automated high-throughput microscope](assets/SMLM_setup_github.png)

# Data analysis

The repository also provides Python scripts for reconstructing objective-scanned 3D SMLM data acquired using the microscope’s 3D imaging beam path.
