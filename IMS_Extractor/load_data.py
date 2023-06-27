import numpy as np
import matplotlib.pyplot as plt
from obs.molecule import Molecule
from obs.id import ID
from obs.envelope import Envelope, Peak

data = np.load("all_data_saved.npy", allow_pickle=True)

num_in_dt = 0
start_rt = data[0]['rt']
for i in range(10000):
	if data[i]['rt'] != start_rt:
		num_in_dt = i
		break

data = data.reshape((int(len(data)/num_in_dt), num_in_dt))
dt_times = [i["dt"] for i in data[0]]
rt_times = [i["rt"] for i in data[:, 0]]

mz_values_m_1 = np.vectorize((lambda x: x["_peaks"][1].mz))(data)
mz_values_m0 = np.vectorize((lambda x: x["_peaks"][1].mz))(data)
mz_values_m1 = np.vectorize((lambda x: x["_peaks"][2].mz))(data)
mz_values_m2 = np.vectorize((lambda x: x["_peaks"][3].mz))(data)
mz_values_mp1 = np.vectorize((lambda x: x["_peaks"][1].mz))(data)

ab_values_m_1 = np.vectorize((lambda x: x["_peaks"][0].ab))(data)
ab_values_m0 = np.vectorize((lambda x: x["_peaks"][1].ab))(data)
ab_values_m1 = np.vectorize((lambda x: x["_peaks"][2].ab))(data)
ab_values_m2 = np.vectorize((lambda x: x["_peaks"][3].ab))(data)
ab_values_mp1 = np.vectorize((lambda x: x["_peaks"][4].ab))(data)

X, Y, = np.meshgrid(dt_times, rt_times)

# Load in the dt ranges into Molecule:
dt_molecule = Molecule("Test Molecule-DT")
all_abunds = np.moveaxis(np.dstack((ab_values_m0, ab_values_m1, ab_values_m2)), -1, 0)
for i in range(all_abunds.shape[1]):
	dt_chrom = all_abunds[:, i, :]
	rt = rt_times[i]
	dt = dt_times
	dt_molecule.add_molecule_id(rt, dt, dt_chrom, 'dt')
dt_molecule.analyze_peaks(5)
dt_molecule.choose_peak()

# Load in the rt ranges into Molecule:
rt_molecule = Molecule("Test Molecule-RT")
for i in range(all_abunds.shape[2]):
	rt_chrom = all_abunds[:, :, i]
	rt = rt_times
	dt = dt_times[i]
	rt_molecule.add_molecule_id(dt, rt, rt_chrom, 'rt')
rt_molecule.analyze_peaks(10)
rt_molecule.choose_peak()
print('hi')

dt_range = dt_molecule.t_range
rt_range = rt_molecule.t_range
ax = plt.axes(projection="3d")
print(dt_times[0])
# ax.plot_surface(X[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
# 				Y[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
# 				ab_values_m0[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
# 				rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.plot_wireframe(X[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
				  Y[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
				  ab_values_m0[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
				  color='black')
ax.plot_wireframe(X[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
				  Y[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
				  ab_values_m1[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
				  color='red')
ax.plot_wireframe(X[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
				  Y[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
				  ab_values_m2[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
				  color='blue')

ax.set_xlabel('DT')
ax.set_ylabel("RT")
ax.set_zlabel("Intensity")
plt.show()


chosen_data = np.stack([ab_values_m_1[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
						ab_values_m0[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
						ab_values_m1[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
						ab_values_m2[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
						ab_values_mp1[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],],
					   0)

mz_data = np.stack([mz_values_m_1[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
					mz_values_m0[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
					mz_values_m1[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
					mz_values_m2[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],
					mz_values_mp1[rt_range[0]: rt_range[1], dt_range[0]: dt_range[1]],],
				   0)
used_rt = rt_times[rt_range[0]: rt_range[1]]
ab_dt_collapsed = chosen_data.sum(2)
mz_dt_collapsed = mz_data.mean(2)
plt.plot(ab_dt_collapsed[1])
plt.plot(ab_dt_collapsed[2])
plt.plot(ab_dt_collapsed[3])
plt.show()

id = ID(rt=np.median(rt_times),
		mz=mz_values_m0.mean(),
		mass=mz_values_m0.mean(),
		z=1,
		n_isos=3)
for scan in range(ab_dt_collapsed.shape[1]):
	envelope_abs = ab_dt_collapsed[:, scan]
	envelope_mzs = mz_dt_collapsed[:, scan]
	rt = used_rt[scan]
	env = Envelope(rt=rt,
				   n_lookback=1,
				   n_lookahead=1)
	for i in range(ab_dt_collapsed.shape[0]):
		p = Peak(envelope_mzs[i],
				 envelope_abs[i],
				 i-1)
		env.append_peak(p)
	id.append_envelope(env)
id.aggregate_envelopes()

print('hi')
