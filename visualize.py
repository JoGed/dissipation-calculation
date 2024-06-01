import os
import numpy as np
from matplotlib import pyplot as plt
from constants import UNIT_NM

def plt_refractive_index(directory, lambd_ccpr, eps_ccpr, timeconvention=1, plot_eps=False, lambd_data=None, eps_data=None, singlePic=True):
  """
  Plot the refractive index or permittivity based on given permittivity within a wavelenght range as an input
  """
  n_ccpr = np.sqrt(eps_ccpr)
  n_data = None
  if (timeconvention == -1): n_ccpr = np.sqrt(np.conjugate(eps_ccpr))
  if (isinstance(eps_data, np.ndarray)): n_data = np.sqrt(eps_data)
  # ---------------------------------------------
  label_data_real = r'$\varepsilon^\prime_{\mathrm{Exp}}$' if plot_eps else r'$n^\prime_{\mathrm{Exp}}$'
  label_data_imag = r'$\varepsilon^{\prime \prime}_{\mathrm{Exp}}$' if plot_eps else r'$n^{\prime \prime}_{\mathrm{Exp}}$'
  label_ccpr_real = r'$\varepsilon^\prime_{\mathrm{CCPR}}$' if plot_eps else r'$n^\prime_{\mathrm{CCPR}}$'
  label_ccpr_imag = r'$\varepsilon^{\prime \prime}_{\mathrm{CCPR}}$' if plot_eps else r'$n^{\prime \prime}_{\mathrm{CCPR}}$'
   # ---------------------------------------------
  
  if singlePic:
    fig, ax = plt.subplots(1, 1, figsize=(3.8, 2.5))
    ax.grid(True, zorder=0, linestyle='dotted')
    if (isinstance(lambd_data, np.ndarray) and isinstance(n_data, np.ndarray)):
        arr_data = n_data**2 if plot_eps else n_data
        ax.scatter(UNIT_NM * lambd_data, np.real(arr_data), label=label_data_real, marker = 'o',facecolor = 'none', edgecolor = 'r', zorder=3)
        ax.scatter(UNIT_NM * lambd_data, np.imag(arr_data), label=label_data_imag, marker = 'o',facecolor = 'none', edgecolor = 'g', zorder=3) 
    arr_ccpr = n_ccpr**2 if plot_eps else n_ccpr
    ax.plot(UNIT_NM * lambd_ccpr, np.real(arr_ccpr[0]), label=label_ccpr_real , color='r')
    ax.plot(UNIT_NM * lambd_ccpr, np.imag(arr_ccpr[0]), label=label_ccpr_imag, color='g')
    ax.set_ylabel(r'$\varepsilon$') if plot_eps else ax.set_ylabel(r'$n$')
    ax.legend()
    ax.set_xlim(UNIT_NM *  lambd_ccpr[0], UNIT_NM *  lambd_ccpr[-1])
    ax.set_xlabel(r'$\lambda\,[\mathrm{nm}]$')
  else:
    fig, axs = plt.subplots(1, 1, figsize=(3.8, 2.5))
    ax.grid(True, zorder=0, linestyle='dotted')
    titles = [r'$\varepsilon_{\mathrm{xx}}$', r'$\varepsilon_{\mathrm{yy}}$', r'$\varepsilon_{\mathrm{zz}}$'] if plot_eps else [r'$n_{\mathrm{xx}}$', r'$n_{\mathrm{yy}}$', r'$n_{\mathrm{zz}}$']
    arr_ccpr = n_ccpr**2 if plot_eps else n_ccpr
    for i, ax in enumerate(axs):
        if (isinstance(lambd_data, np.ndarray) and isinstance(n_data, np.ndarray)):
          arr_data = n_data**2 if plot_eps else n_data
          ax.scatter(UNIT_NM * lambd_data, np.real(arr_data), label=label_data_real, marker = 'o',facecolor = 'none', edgecolor = 'r', zorder=3)
          ax.scatter(UNIT_NM * lambd_data, np.imag(arr_data), label=label_data_imag, marker = 'o',facecolor = 'none', edgecolor = 'g', zorder=3)
        ax.plot(UNIT_NM * lambd_ccpr, np.real(arr_ccpr[i]), label=label_ccpr_real, color='r')
        ax.plot(UNIT_NM * lambd_ccpr, np.imag(arr_ccpr[i]), label=label_ccpr_imag, color='g')
        ax.set_title(titles[i])
        ax.set_ylabel(r'$\varepsilon$') if plot_eps else ax.set_ylabel(r'$n$')
        ax.legend()
      
        ax.set_xlim(UNIT_NM *  lambd_ccpr[0], UNIT_NM *  lambd_ccpr[-1])
        ax.set_xlabel(r'$\lambda\,[\mathrm{nm}]$')
   # ---------------------------------------------

  file_name = "eps.pdf" if plot_eps else "refractive_index.pdf"
  plt.legend(loc='upper left', bbox_to_anchor=(1,1), framealpha=1)
  plt.tight_layout()
  plt.savefig( os.path.join(directory, file_name), dpi = 800)
  plt.show()

  return 1


def plt_qcoeffs(directory, lambd, qext_mie, qsca_mie, qabs_mie, lambd_data, qext_data, qsca_data, qabs_data, lambd_integr, qabs_integr):
  """
  Plot the Q coefficients within a given wavelenght range
  """
  plot_data = [
    (r'$Q_{\mathrm{sca}}$', 'blue', lambd, qsca_mie, lambd_data, qsca_data),
    (r'$Q_{\mathrm{abs}}$', 'red', lambd, qabs_mie, lambd_data, qabs_data),
  ]
  fig, ax = plt.subplots(1, 1, figsize=(6.5, 3.5))
  ax.grid(True, zorder=0, linestyle='dotted')
  for label, color, lambda_mie, q_mie, lambda_data, q_data in plot_data:
      ax.plot(UNIT_NM * lambda_mie, q_mie, label=f'{label}' + r'$_\mathrm{,\,Mie, analytical}$', color=color, linestyle='dotted')
      ax.plot(UNIT_NM * lambda_data, q_data, label=f'{label}' + r'$_\mathrm{,\,FDTD,\,Surface.integr.}$', color=color)
  ax.plot(UNIT_NM * lambd_integr, qabs_integr, label=r'$Q_{\mathrm{abs,\, FDTD,\, Vol.integr.}}$', color="g")
  ax.set_xlabel(r'$\lambda\,[\mathrm{nm}]$')
  ax.set_ylabel(r'$Q$')
  ax.legend()
  ax.set_xlim(UNIT_NM * lambd[0], UNIT_NM * lambd[-1])
  ax.set_ylim(0,)
  plt.tight_layout()
  plt.savefig(os.path.join(directory, 'Qs.pdf'))
  plt.show()

  return 1

def plt_quantity(directory, lambd, quantity, q_str):
  """
  Plot a quantity within a given wavelenght range
  """
  lambd_nm = UNIT_NM * lambd
  fig, ax = plt.subplots(1, 1, figsize=(3.8, 3.2))
  ax.grid(True, zorder=0, linestyle='dotted')
  ax.plot(lambd_nm, quantity)
  ax.set_ylabel(q_str)
  ax.set_xlabel(r'$\lambda\,[\mathrm{nm}]$')
  ax.set_xlim(lambd_nm[0], lambd_nm[-1])
  plt.tight_layout()
  plt.savefig(os.path.join(directory, q_str + '.pdf'))
  plt.show()
  return 1

