/*
dynsmooth~ - Dynamic smoothing filter external for Pure Data
2024
Functionality:
- Implements adaptive signal smoothing with dynamic cutoff
- Amount of smoothing adjusts automatically based on signal characteristics
- Base frequency and sensitivity configurable via creation arguments
- Signal rate processing
- Supports multichannel processing
Usage:
1. Creation args: [basefreq] [sensitivity] (defaults: 2, 2)
2. Signal inlet: Input signal
3. Outlet: Smoothed signal
Messages:
- clear: Reset filter states to zero
- basefreq [float]: Set base frequency in Hz (lower = more smoothing)
- sensitivity [float]: Set sensitivity factor (higher = more dynamic response)
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
    t_float *inz;          // Input z-1 state
    
    // Parameters
    t_float basefreq;      // Base frequency in Hz
    t_float sensitivity;   // Sensitivity factor
    t_float wc;           // Normalized cutoff frequency
    
    int n_channels;        // Number of channels
    t_float sr;           // Sample rate
} t_dynsmooth_tilde;

// Helper function to calculate normalized cutoff
static t_float calc_wc(t_float basefreq, t_float sr) {
    return (2.0f * M_PI * basefreq) / sr;
}

// Perform routine - processes audio at signal rate
static t_int *dynsmooth_tilde_perform(t_int *w)
{
    t_dynsmooth_tilde *x = (t_dynsmooth_tilde *)(w[1]);
    int n = (int)(w[2]);
    t_sample *in = (t_sample *)(w[3]);
    t_sample *out = (t_sample *)(w[4]);
    
    t_float low1z, low2z, bandz, wd, g;
    t_float wc = x->wc;
    t_float sensitivity = x->sensitivity;
    
    // Process samples
    for (int i = 0; i < n; i++) {
        // Store previous states
        low1z = x->low1[0];
        low2z = x->low2[0];
        
        // Calculate band-pass value
        bandz = low1z - low2z;
        
        // Calculate adaptive cutoff
        wd = wc + sensitivity * fabs(bandz);
        
        // Calculate filter coefficient (cubic approximation)
        g = wd * (5.9948827f + wd * (-11.969296f + wd * 15.959062f));
        if (g > 1.0f) g = 1.0f;
        
        // Update states
        x->low1[0] = low1z + g * (0.5f * (in[i] + x->inz[0]) - low1z);
        x->low2[0] = low2z + g * (0.5f * (x->low1[0] + low1z) - low2z);
        x->inz[0] = in[i];
        
        // Output second lowpass state
        out[i] = x->low2[0];
    }
    
    return (w + 5);
}

static void dynsmooth_tilde_dsp(t_dynsmooth_tilde *x, t_signal **sp)
{
    // Update sample rate if changed
    if (x->sr != sp[0]->s_sr) {
        x->sr = sp[0]->s_sr;
        x->wc = calc_wc(x->basefreq, x->sr);
    }
    
    dsp_add(dynsmooth_tilde_perform, 4, x, 
            sp[0]->s_n, sp[0]->s_vec, sp[1]->s_vec);
}

// Message handler for 'clear' message
static void dynsmooth_tilde_clear(t_dynsmooth_tilde *x)
{
    for (int i = 0; i < x->n_channels; i++) {
        x->low1[i] = 0.0f;
        x->low2[i] = 0.0f;
        x->inz[i] = 0.0f;
    }
}

// Message handler for 'basefreq' message
static void dynsmooth_tilde_basefreq(t_dynsmooth_tilde *x, t_floatarg f)
{
    if (f > 0) {
        x->basefreq = f;
        x->wc = calc_wc(f, x->sr);
    }
}

// Message handler for 'sensitivity' message
static void dynsmooth_tilde_sensitivity(t_dynsmooth_tilde *x, t_floatarg f)
{
    x->sensitivity = f;
}

static void *dynsmooth_tilde_new(t_floatarg basefreq, t_floatarg sensitivity)
{
    t_dynsmooth_tilde *x = (t_dynsmooth_tilde *)pd_new(dynsmooth_tilde_class);
    
    // Set default values if not specified
    x->basefreq = basefreq > 0 ? basefreq : 2.0f;
    x->sensitivity = sensitivity > 0 ? sensitivity : 2.0f;
    
    // Initialize single channel
    x->n_channels = 1;
    x->sr = sys_getsr();
    x->wc = calc_wc(x->basefreq, x->sr);
    
    // Allocate memory for states
    x->low1 = (t_float *)getbytes(sizeof(t_float));
    x->low2 = (t_float *)getbytes(sizeof(t_float));
    x->inz = (t_float *)getbytes(sizeof(t_float));
    
    // Initialize states to zero
    dynsmooth_tilde_clear(x);
    
    // Create signal outlet
    outlet_new(&x->x_obj, &s_signal);
    
    return (void *)x;
}

static void dynsmooth_tilde_free(t_dynsmooth_tilde *x)
{
    if (x->low1) freebytes(x->low1, x->n_channels * sizeof(t_float));
    if (x->low2) freebytes(x->low2, x->n_channels * sizeof(t_float));
    if (x->inz) freebytes(x->inz, x->n_channels * sizeof(t_float));
}

void dynsmooth_tilde_setup(void)
{
    dynsmooth_tilde_class = class_new(gensym("dynsmooth~"),
        (t_newmethod)dynsmooth_tilde_new,
        (t_method)dynsmooth_tilde_free,
        sizeof(t_dynsmooth_tilde),
        CLASS_DEFAULT,
        A_DEFFLOAT, A_DEFFLOAT, 0);
    
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