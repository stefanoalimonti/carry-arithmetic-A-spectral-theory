/*
 * A14_isochrone_mc.c
 *
 * Monte Carlo computation of the isochrone profile f(s) = <c_{M-1} - 1>(s)
 * where s = arctan(X) + arctan(Y) in the D-odd angular triangle.
 *
 * Uses same carry chain logic as f81: full carries array, M_top finding,
 * ULC check, cm1 = carries[M_top - 1].
 *
 * Usage: ./A14_isochrone_mc [d] [n_samples] [n_bins]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

static uint64_t s_rng[4];
static uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}
static uint64_t next_rand(void) {
    const uint64_t result = rotl(s_rng[1] * 5, 7) * 9;
    const uint64_t t = s_rng[1] << 17;
    s_rng[2] ^= s_rng[0]; s_rng[3] ^= s_rng[1];
    s_rng[1] ^= s_rng[2]; s_rng[0] ^= s_rng[3];
    s_rng[2] ^= t; s_rng[3] = rotl(s_rng[3], 45);
    return result;
}
static void seed_rng(uint64_t seed) {
    s_rng[0] = seed; s_rng[1] = seed ^ 0x9E3779B97F4A7C15ULL;
    s_rng[2] = seed ^ 0x6C62272E07BB0142ULL;
    s_rng[3] = seed ^ 0xBF58476D1CE4E5B9ULL;
    for (int i = 0; i < 20; i++) next_rand();
}

int main(int argc, char **argv) {
    int d = 64;
    long long n_samples = 200000000LL;
    int n_bins = 50;

    if (argc > 1) d = atoi(argv[1]);
    if (argc > 2) n_samples = atoll(argv[2]);
    if (argc > 3) n_bins = atoi(argv[3]);

    seed_rng(0xCAFEBABE93ULL);

    int D_max = 2 * d;
    int max_pos = D_max + 20;
    double pi4 = M_PI / 4.0;
    double bin_width = pi4 / n_bins;

    double *sum_cm1 = calloc(n_bins, sizeof(double));
    double *sum_cm1_sq = calloc(n_bins, sizeof(double));
    long long *count = calloc(n_bins, sizeof(long long));

    double *sum_cm1_even = calloc(n_bins, sizeof(double));
    long long *count_even = calloc(n_bins, sizeof(long long));

    long long n_ulc = 0, n_even = 0, n_odd = 0;
    double total_c1_num = 0;
    long long total_c1_den = 0;

    printf("A14: ISOCHRONE PROFILE f(s) = <c_{M-1} - 1>(s)\n");
    printf("d = %d, n_samples = %lld, n_bins = %d\n", d, n_samples, n_bins);
    printf("bin_width = %.6f (pi/4 / %d)\n\n", bin_width, n_bins);

    int *carries = malloc((max_pos + 10) * sizeof(int));

    for (long long iter = 0; iter < n_samples; iter++) {
        if (iter > 0 && iter % 50000000LL == 0) {
            double c1_tmp = (total_c1_den > 0) ? total_c1_num / total_c1_den - 1.0 : 0;
            printf("  Progress: %lld / %lldM, ULC=%lld, c1~%.8f\n",
                   iter / 1000000, n_samples / 1000000, n_ulc, c1_tmp);
            fflush(stdout);
        }

        uint64_t p_val = (1ULL << (d-1)) | (next_rand() & ((1ULL << (d-1)) - 1));
        uint64_t q_val = (1ULL << (d-1)) | (next_rand() & ((1ULL << (d-1)) - 1));

        memset(carries, 0, (max_pos + 10) * sizeof(int));
        int last_product_bit = 0;

        for (int pos = 0; pos < D_max; pos++) {
            int conv = 0;
            int lo = (pos >= d) ? pos - d + 1 : 0;
            int hi = (pos < d) ? pos : d - 1;
            for (int i = lo; i <= hi; i++) {
                int gi = (int)((p_val >> i) & 1);
                int hj = (int)((q_val >> (pos - i)) & 1);
                conv += gi * hj;
            }
            int total = conv + carries[pos];
            int bit = total & 1;
            carries[pos + 1] = total >> 1;
            if (bit) last_product_bit = pos;
        }
        for (int pos = D_max; pos < max_pos; pos++) {
            if (carries[pos] == 0) break;
            int total = carries[pos];
            int bit = total & 1;
            carries[pos + 1] = total >> 1;
            if (bit) last_product_bit = pos;
        }

        int D_bits = last_product_bit + 1;

        int M_top = 0;
        for (int pos = max_pos; pos >= 1; pos--) {
            if (carries[pos] > 0) { M_top = pos; break; }
        }

        if (M_top <= 0 || carries[M_top] != 1) continue;
        n_ulc++;

        int cm1 = carries[M_top - 1];
        total_c1_num += cm1;
        total_c1_den++;

        double X = (double)(p_val - (1ULL << (d-1))) / (double)(1ULL << (d-1));
        double Y = (double)(q_val - (1ULL << (d-1))) / (double)(1ULL << (d-1));
        double s_val = atan(X) + atan(Y);

        if (D_bits % 2 != 0) {
            n_odd++;
            int bin = (int)(s_val / bin_width);
            if (bin >= n_bins) bin = n_bins - 1;
            if (bin < 0) bin = 0;

            double val = cm1 - 1.0;
            sum_cm1[bin] += val;
            sum_cm1_sq[bin] += val * val;
            count[bin]++;
        } else {
            n_even++;
            double s_even = s_val - pi4;
            if (s_even < 0) s_even = 0;
            int bin = (int)(s_even / bin_width);
            if (bin >= n_bins) bin = n_bins - 1;
            sum_cm1_even[bin] += cm1 - 1.0;
            count_even[bin]++;
        }
    }

    double c1 = total_c1_num / total_c1_den - 1.0;
    printf("\n{'='*70}\n");
    printf("=== GLOBAL RESULTS ===\n");
    printf("Total samples: %lld\n", n_samples);
    printf("ULC valid: %lld (%.4f%%)\n", n_ulc, 100.0 * n_ulc / n_samples);
    printf("D-even: %lld (%.6f), D-odd: %lld (%.6f)\n",
           n_even, (double)n_even / n_ulc, n_odd, (double)n_odd / n_ulc);
    printf("c1 = E[cm1]-1 = %.10f\n", c1);
    printf("pi/18         = %.10f\n", M_PI / 18.0);
    printf("difference    = %.2e\n\n", c1 - M_PI / 18.0);

    printf("=== D-ODD ISOCHRONE PROFILE ===\n");
    printf("%-6s %-10s %-10s %-14s %-14s %-12s %-10s\n",
           "bin", "s_mid", "T=tan(s)", "f(s)=<cm1-1>", "f_err", "count", "I(s)");
    printf("------------------------------------------------------------------------\n");

    double integral_fI = 0.0;
    double integral_I = 0.0;

    for (int b = 0; b < n_bins; b++) {
        double s_mid = (b + 0.5) * bin_width;
        double T = tan(s_mid);

        double I_s;
        if (T < 1e-10) {
            I_s = 2.0 / 3.0;
        } else {
            I_s = 2.0 * (1 + T * T) * (T * T - log(1 + T * T)) / (T * T * T);
        }

        double f_s = 0.0, f_err = 0.0;
        if (count[b] > 0) {
            f_s = sum_cm1[b] / count[b];
            double var = sum_cm1_sq[b] / count[b] - f_s * f_s;
            if (var > 0) f_err = sqrt(var / count[b]);
        }

        integral_I += I_s * bin_width;
        integral_fI += f_s * I_s * bin_width;

        printf("%-6d %-10.6f %-10.6f %-14.8f %-14.8f %-12lld %-10.6f\n",
               b, s_mid, T, f_s, f_err, count[b], I_s);
    }

    printf("\n=== INTEGRAL CHECKS ===\n");
    printf("∫ I(s) ds [0,pi/4]     = %.10f  (theory P_o = %.10f)\n",
           integral_I, 2 * log(2) - 1);
    printf("∫ f(s)*I(s) ds [0,pi/4]= %.10f  (= Sigma_odd)\n", integral_fI);

    double Sigma_even = 1.0 + 3 * log(0.75);
    double c1_from_integral = integral_fI + Sigma_even;
    printf("c1 = Sigma_even + Sigma_odd = %.10f + %.10f = %.10f\n",
           Sigma_even, integral_fI, c1_from_integral);

    double target_odd = M_PI / 18.0 - Sigma_even;
    printf("Target Sigma_odd = %.10f\n", target_odd);
    printf("pi/18            = %.10f\n", M_PI / 18.0);

    printf("\n=== PROFILE ANALYSIS ===\n");

    double sum_f_w = 0, sum_sf_w = 0, sum_s2_w = 0, sum_s_w = 0;
    long long total_w = 0;
    for (int b = 0; b < n_bins; b++) {
        if (count[b] < 100) continue;
        double s_mid = (b + 0.5) * bin_width;
        double f_s = sum_cm1[b] / count[b];
        double w = (double)count[b];
        sum_f_w += w * f_s;
        sum_sf_w += w * s_mid * f_s;
        sum_s2_w += w * s_mid * s_mid;
        sum_s_w += w * s_mid;
        total_w += count[b];
    }

    if (total_w > 0) {
        double mean_f = sum_f_w / total_w;
        double mean_s = sum_s_w / total_w;
        double cov_sf = sum_sf_w / total_w - mean_s * mean_f;
        double var_s = sum_s2_w / total_w - mean_s * mean_s;
        double slope = cov_sf / var_s;
        double intercept = mean_f - slope * mean_s;

        printf("Linear fit f(s) ≈ a + b*s:\n");
        printf("  a = %.8f, b = %.8f\n", intercept, slope);
        printf("  f(0) ~ %.8f, f(pi/4) ~ %.8f\n",
               intercept, intercept + slope * pi4);
    }

    printf("\n=== ANALYTICAL FIT TESTS ===\n");
    printf("Testing f(s) = C*s, f(s) = C*arctan(tan(s)), f(s) = C*tan(s):\n\n");
    for (int b = 0; b < n_bins; b += 5) {
        if (count[b] < 100) continue;
        double s_mid = (b + 0.5) * bin_width;
        double T = tan(s_mid);
        double f_s = sum_cm1[b] / count[b];
        printf("  s=%.4f T=%.4f: f=%.8f  f/s=%.6f  f/T=%.6f  f/atan(T)=%.6f  f/T^2=%.6f\n",
               s_mid, T, f_s, f_s/s_mid, (T>1e-10)?f_s/T:0,
               (s_mid>1e-10)?f_s/s_mid:0, (T>1e-10)?f_s/(T*T):0);
    }

    free(sum_cm1); free(sum_cm1_sq); free(count);
    free(sum_cm1_even); free(count_even); free(carries);
    return 0;
}
