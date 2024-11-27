/*
dynsmooth~ - Dynamic smoothing filter external for Pure Data
Ben Wesch, 2024
Functionality:
- Implements adaptive signal smoothing with dynamic cutoff
- Two modes of operation: accurate (default) and efficient
- Base frequency and sensitivity configurable via creation arguments
- Signal rate processing
Usage:
1. Creation args: [-e] [basefreq] [sensitivity]
   -e: Optional flag for efficient mode (default is accurate mode)
   basefreq: Base frequency in Hz (default: 2)
   sensitivity: Sensitivity factor (default: accurate=2, efficient=0.5)
2. Signal inlet: Input signal
3. Outlet: Smoothed signal
Messages:
- clear: Reset filter states to zero
- basefreq [float]: Set base frequency in Hz (lower = more smoothing)
- sensitivity [float]: Set sensitivity factor
*/

#include "m_pd.h"
#include <math.h>

static t_class *dynsmooth_tilde_class;

typedef struct _dynsmooth_tilde {
    t_object x_obj;
    t_sample f_dummy;      // Dummy float for main signal inlet
    
    // Filter state per channel
    t_float *low1;         // First lowpass state
    t_float *low2;         // Second lowpass state
    t_float *inz;          // Input z-1 state (only used in accurate mode)
    
    // Parameters
    t_float basefreq;      // Base frequency in Hz
    t_float sensitivity;   // Sensitivity factor
    t_float wc;           // Normalized cutoff frequency
    t_float gc;           // Tan-based coefficient (efficient mode)
    t_float g0;           // Base coefficient (efficient mode)
    t_float sense;        // Scaled sensitivity (efficient mode)
    
    int n_channels;        // Number of channels
    t_float sr;           // Sample rate
    int efficient_mode;    // Mode flag (0 = accurate, 1 = efficient)
} t_dynsmooth_tilde;

// Helper functions
static t_float calc_wc(t_float basefreq, t_float sr) {
    return (2.0f * M_PI * basefreq) / sr;
}

static void update_efficient_coeffs(t_dynsmooth_tilde *x) {
    x->gc = tanf(M_PI * x->wc);
    x->g0 = 2.0f * x->gc / (1.0f + x->gc);
    x->sense = x->sensitivity * 4.0f;
}

// Perform routine for accurate mode
static t_int *dynsmooth_tilde_perform_accurate(t_int *w)
{
    t_dynsmooth_tilde *x = (t_dynsmooth_tilde *)(w[1]);
    int n = (int)(w[2]);
    t_sample *in = (t_sample *)(w[3]);
    t_sample *out = (t_sample *)(w[4]);
    
    t_float low1z, low2z, bandz, wd, g;
    t_float wc = x->wc;
    t_float sensitivity = x->sensitivity;

    for (int c = 0; c < x->n_channels; c++) {
        for (int i = 0; i < n; i++) {
            // Store previous states
            low1z = x->low1[c];
            low2z = x->low2[c];

            // Calculate band-pass value
            bandz = low1z - low2z;
            
            // Calculate adaptive cutoff
            wd = wc + sensitivity * fabs(bandz);

            // Calculate filter coefficient (cubic approximation)
            g = wd * (5.9948827f + wd * (-11.969296f + wd * 15.959062f));
            if (g > 1.0f) g = 1.0f;
            
            // Update states
            x->low1[c] = low1z + g * (0.5f * (in[i + c * n] + x->inz[c]) - low1z);
            x->low2[c] = low2z + g * (0.5f * (x->low1[c] + low1z) - low2z);
            x->inz[c] = in[i + c * n];
            
            // Output second lowpass state
            out[i + c * n] = x->low2[c];
        }
    }

    return (w + 5);
}

// Perform routine for efficient mode
static t_int *dynsmooth_tilde_perform_efficient(t_int *w)
{
    t_dynsmooth_tilde *x = (t_dynsmooth_tilde *)(w[1]);
    int n = (int)(w[2]);
    t_sample *in = (t_sample *)(w[3]);
    t_sample *out = (t_sample *)(w[4]);
    
    t_float low1z, low2z, bandz, g;
    t_float g0 = x->g0;
    t_float sense = x->sense;
    
    for (int c = 0; c < x->n_channels; c++) {
        for (int i = 0; i < n; i++) {
            low1z = x->low1[c];
            low2z = x->low2[c];
            bandz = low1z - low2z;
            
            g = g0 + sense * fabs(bandz);
            if (g > 1.0f) g = 1.0f;
            
            x->low1[c] = low1z + g * (in[i + c * n] - low1z);
            x->low2[c] = low2z + g * (x->low1[c] - low2z);
            
            out[i + c * n] = x->low2[c];
        }
    }    
    return (w + 5);
}

static void dynsmooth_tilde_dsp(t_dynsmooth_tilde *x, t_signal **sp)
{
    // Update sample rate if changed
    if (x->sr != sp[0]->s_sr) {
        x->sr = sp[0]->s_sr;
        x->wc = calc_wc(x->basefreq, x->sr);
        if (x->efficient_mode) {
            update_efficient_coeffs(x);
        }
    }

    int n_channels = sp[0]->s_nchans;
    // Reallocate memory
    x->low1 = (t_float *)resizebytes(x->low1, sizeof(t_float) * x->n_channels, sizeof(t_float) * n_channels);
    x->low2 = (t_float *)resizebytes(x->low2, sizeof(t_float) * x->n_channels, sizeof(t_float) * n_channels);
    if (!x->efficient_mode) {
        x->inz = (t_float *)resizebytes(x->inz, sizeof(t_float) * x->n_channels, sizeof(t_float) * n_channels);
    }

    for (int i = x->n_channels; i < n_channels; i++) {
        x->low1[i] = 0.0f;
        x->low2[i] = 0.0f;
        if (!x->efficient_mode) x->inz[i] = 0.0f;
    }

    x->n_channels = n_channels;
    // Set channel count based on input
    signal_setmultiout(&sp[1], n_channels);

    // Select perform routine based on mode
    if (x->efficient_mode) {
        dsp_add(dynsmooth_tilde_perform_efficient, 4, x, 
                sp[0]->s_n, sp[0]->s_vec, sp[1]->s_vec);
    } else {
        dsp_add(dynsmooth_tilde_perform_accurate, 4, x, 
                sp[0]->s_n, sp[0]->s_vec, sp[1]->s_vec);
    }
}

static void dynsmooth_tilde_clear(t_dynsmooth_tilde *x)
{
    for (int i = 0; i < x->n_channels; i++) {
        x->low1[i] = 0.0f;
        x->low2[i] = 0.0f;
        if (!x->efficient_mode) {
            x->inz[i] = 0.0f;
        }
    }
}

static void dynsmooth_tilde_basefreq(t_dynsmooth_tilde *x, t_floatarg f)
{
    if (f > 0) {
        x->basefreq = f;
        x->wc = calc_wc(f, x->sr);
        if (x->efficient_mode) {
            update_efficient_coeffs(x);
        }
    }
}

static void dynsmooth_tilde_sensitivity(t_dynsmooth_tilde *x, t_floatarg f)
{
    x->sensitivity = f;
    if (x->efficient_mode) {
        x->sense = f * 4.0f;
    }
}

static void *dynsmooth_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    (void)s;

    t_dynsmooth_tilde *x = (t_dynsmooth_tilde *)pd_new(dynsmooth_tilde_class);
    
    // Parse arguments
    x->efficient_mode = 0;  // Default to accurate mode
    x->basefreq = 2.0f;    // Default base frequency
    x->sensitivity = 2.0f;  // Default sensitivity for accurate mode
    
    // Check for efficient mode flag
    if (argc && argv->a_type == A_SYMBOL && atom_getsymbol(argv) == gensym("-e")) {
        x->efficient_mode = 1;
        x->sensitivity = 0.5f; // Different default for efficient mode
        argc--, argv++;
    }
    
    // Get base frequency if specified
    if (argc) {
        float basefreq = atom_getfloat(argv);
        if (basefreq > 0) x->basefreq = basefreq;
        argc--, argv++;
    }

    // Get sensitivity if specified
    if (argc) {
        float sens = atom_getfloat(argv);
        if (sens > 0) x->sensitivity = sens;
        argc--, argv++;
    }
    
    // Initialize
    x->n_channels = 1;
    x->sr = sys_getsr();
    x->wc = calc_wc(x->basefreq, x->sr);
    
    if (x->efficient_mode) {
        update_efficient_coeffs(x);
    }
    
    // Allocate memory
    x->low1 = (t_float *)getbytes(sizeof(t_float));
    x->low2 = (t_float *)getbytes(sizeof(t_float));
    if (!x->efficient_mode) {
        x->inz = (t_float *)getbytes(sizeof(t_float));
    }
    
    // Initialize states
    dynsmooth_tilde_clear(x);
    
    // Create signal outlet
    outlet_new(&x->x_obj, &s_signal);
    
    return (void *)x;
}

static void dynsmooth_tilde_free(t_dynsmooth_tilde *x)
{
    if (x->low1) freebytes(x->low1, x->n_channels * sizeof(t_float));
    if (x->low2) freebytes(x->low2, x->n_channels * sizeof(t_float));
    if (!x->efficient_mode && x->inz) {
        freebytes(x->inz, x->n_channels * sizeof(t_float));
    }
}

void dynsmooth_tilde_setup(void)
{
    dynsmooth_tilde_class = class_new(gensym("dynsmooth~"),
        (t_newmethod)dynsmooth_tilde_new,
        (t_method)dynsmooth_tilde_free,
        sizeof(t_dynsmooth_tilde),
        CLASS_MULTICHANNEL,
        A_GIMME, 0);  // Changed to A_GIMME to handle optional flag
    
    class_addmethod(dynsmooth_tilde_class,
        (t_method)dynsmooth_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(dynsmooth_tilde_class,
        (t_method)dynsmooth_tilde_clear, gensym("clear"), 0);
    class_addmethod(dynsmooth_tilde_class,
        (t_method)dynsmooth_tilde_basefreq, gensym("basefreq"), A_FLOAT, 0);
    class_addmethod(dynsmooth_tilde_class,
        (t_method)dynsmooth_tilde_sensitivity, gensym("sensitivity"), A_FLOAT, 0);
    
    CLASS_MAINSIGNALIN(dynsmooth_tilde_class, t_dynsmooth_tilde, f_dummy);
}
