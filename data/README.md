# LIQUID EDGE IEEE 802.11ad PHY Simulator - Data

![LIQUID EDGE Logo](../doc/liquid_edge_logo28.png)

> These programs are a part of the system used for the LIQUID EDGE PRIN 2017 project demonstrator.

This folder holds auxiliary data for the simulator operation.

If you want to use the Brno University of Technology (BUT) channel model [2], you have to copy here the folder named *measured_channels* located in the relevant [GitHub repository](https://github.com/jirimilos/802.11ad-phy-sim) of [3]. The final structure should be like this:

* data
  * measured_channels
    * CIR.mat
    * CIR2.mat
    * CTF.mat
    * CTF2.mat
    * show_meas_data.m

For instructions on the simulator operation, refer to the main [README](../README.md) file.

## References

1. "IEEE Standard for Information technology--Telecommunications and information exchange between systems--Local and metropolitan area networks--Specific requirements-Part 11: Wireless LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications Amendment 3: Enhancements for Very High Throughput in the 60 GHz Band," in _IEEE Std 802.11ad-2012 (Amendment to IEEE Std 802.11-2012, as amended by IEEE Std 802.11ae-2012 and IEEE Std 802.11aa-2012)_, pp.1-628, 28 Dec. 2012, doi: 10.1109/IEEESTD.2012.6392842
2. P. Liu, J. Blumenstein, N. S. Perović, M. Di Renzo and A. Springer, "Performance of Generalized Spatial Modulation MIMO Over Measured 60GHz Indoor Channels," in _IEEE Transactions on Communications_, vol. 66, no. 1, pp. 133-148, Jan. 2018. doi: 10.1109/TCOMM.2017.2754280 URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8046024&isnumber=8258580
3. J. Blumenstein, J. Milos, L. Polak and C. Mecklenbräuker, "IEEE 802.11ad SC-PHY Layer Simulator: Performance in Real-world 60 GHz Indoor Channels," _2019 IEEE Nordic Circuits and Systems Conference (NORCAS): NORCHIP and International Symposium of System-on-Chip (SoC)_, Helsinki, Finland, 2019, pp. 1-4, doi: 10.1109/NORCHIP.2019.8906960
4. G. Baruffa, L. Rugini, "Improved Channel Estimation and Equalization for Single-Carrier IEEE 802.11ad Receivers," submitted to _Radioengineering_, 2023