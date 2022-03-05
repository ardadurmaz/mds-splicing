# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 14:00:42 2021

@author: durmaz
"""




from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
from sklearn.preprocessing import scale
import tensorflow as tf

def cust_loss(y_true, y_pred):
    
    y_pred = tf.clip_by_value(y_pred, tf.keras.backend.epsilon(), 1.0-tf.keras.backend.epsilon())
    log_loss = tf.math.add(tf.math.multiply(y_true, tf.math.log(y_pred)), tf.math.multiply(tf.math.multiply(1.0-y_true, tf.math.log(1.0-y_pred)), 0.1))
    
    return tf.math.reduce_sum(tf.math.reduce_mean(-1.0*log_loss, axis=-1), axis=-1)

def cust_loss_mse(y_true, y_pred):
    mse_err = tf.keras.losses.mean_squared_error(y_true, y_pred)
    
    return tf.math.reduce_sum(mse_err, axis=-1)


class Sampling(tf.keras.layers.Layer):

    def call(self, inputs):
        z_mean, z_log_var = inputs
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.keras.backend.random_normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon
        

class VariationalAutoEncoder(tf.keras.Model):

    def __init__(self, enc, dec, name="autoencoder"):
        super(VariationalAutoEncoder, self).__init__(name=name)
        self.encoder = enc
        self.decoder = dec
        self.custom_weight = tf.Variable(name='kl_reg', initial_value=0.0, trainable=False, dtype="float32")

    def call(self, inputs):
        # Reconstruction Loss
        z_mean, z_log_var, z, lat_expr_sec, sd_expr_sec, z_expr, lat_mut_sec, sd_mut_sec, z_mut, lat_psi_sec, sd_psi_sec, z_psi = self.encoder(inputs)
        reconstructed = self.decoder([z, z_expr, z_mut, z_psi])
        
        # Add KL divergence regularization loss.
        kl_loss = -0.5 * tf.reduce_mean(
            z_log_var - tf.square(z_mean) - tf.exp(z_log_var) + 1.0
        )
        
        kl_loss_expr = -0.5 * tf.reduce_mean(
            sd_expr_sec - tf.square(lat_expr_sec) - tf.exp(sd_expr_sec) + 1.0
        )
        
        kl_loss_mut = -0.5 * tf.reduce_mean(
            sd_mut_sec - tf.square(lat_mut_sec) - tf.exp(sd_mut_sec) + 1.0
        )
        
        kl_loss_psi = -0.5 * tf.reduce_mean(
            sd_psi_sec - tf.square(lat_psi_sec) - tf.exp(sd_psi_sec) + 1.0
        )
        
        kl_loss_comb = tf.math.reduce_sum([kl_loss, kl_loss_expr, kl_loss_mut, kl_loss_psi])
        
        self.add_loss(self.custom_weight*kl_loss_comb)
        self.add_metric(self.custom_weight*kl_loss_comb, name='kl_loss')
        
        return reconstructed


class EncoderModel(tf.keras.Model):

    def __init__(self, latent_dim, name="encoder"):
        super(EncoderModel, self).__init__(name=name)
        self.z_expr = tf.keras.layers.Dense(latent_dim, activation='linear', name='expr_lat')
        self.z_mut = tf.keras.layers.Dense(latent_dim, activation='linear', name='mut_lat')
        self.z_psi = tf.keras.layers.Dense(latent_dim, activation='linear', name='psi_lat')
        self.z_expr_sd = tf.keras.layers.Dense(latent_dim, activation='linear', name='expr_sd')
        self.z_mut_sd = tf.keras.layers.Dense(latent_dim, activation='linear',  name='mut_sd')
        self.z_psi_sd = tf.keras.layers.Dense(latent_dim, activation='linear', name='psi_sd')
        
        self.z_expr_sec = tf.keras.layers.Dense(latent_dim, activation='linear', name='expr_lat')
        self.z_mut_sec = tf.keras.layers.Dense(latent_dim, activation='linear', name='mut_lat')
        self.z_psi_sec = tf.keras.layers.Dense(latent_dim, activation='linear', name='psi_lat')
        self.z_expr_sd_sec = tf.keras.layers.Dense(latent_dim, activation='linear', name='expr_sd')
        self.z_mut_sd_sec = tf.keras.layers.Dense(latent_dim, activation='linear',  name='mut_sd')
        self.z_psi_sd_sec = tf.keras.layers.Dense(latent_dim, activation='linear', name='psi_sd')

        self.sampling = Sampling()

    def call(self, inputs):
        lat_expr = self.z_expr(inputs[0])
        lat_mut = self.z_mut(inputs[1])
        lat_psi = self.z_psi(inputs[2])
        
        sd_expr = self.z_expr_sd(inputs[0])
        sd_mut = self.z_mut_sd(inputs[1])
        sd_psi = self.z_psi_sd(inputs[2])

        lat_expr_sec = self.z_expr_sec(inputs[0])
        lat_mut_sec = self.z_mut_sec(inputs[1])
        lat_psi_sec = self.z_psi_sec(inputs[2])

        sd_expr_sec = self.z_expr_sd_sec(inputs[0])
        sd_mut_sec = self.z_mut_sd_sec(inputs[1])
        sd_psi_sec = self.z_psi_sd_sec(inputs[2])

        
        z_mean = tf.math.add(tf.math.add(lat_expr, lat_mut), lat_psi)
        z_log_var = tf.math.add(tf.math.add(sd_expr, sd_mut), sd_psi)
        z = self.sampling([z_mean, z_log_var])
        z_expr = self.sampling([lat_expr_sec, sd_expr_sec])
        z_mut = self.sampling([lat_mut_sec, sd_mut_sec])
        z_psi = self.sampling([lat_psi_sec, sd_psi_sec])
        
        return z_mean, z_log_var, z, lat_expr_sec, sd_expr_sec, z_expr, lat_mut_sec, sd_mut_sec, z_mut, lat_psi_sec, sd_psi_sec, z_psi

class DecoderModel(tf.keras.Model):

    def __init__(self, EXPR_DIM, MUT_DIM, PSI_DIM, name="decoder"):
        super(DecoderModel, self).__init__(name=name)
        self.out_expr = tf.keras.layers.Dense(EXPR_DIM, activation='linear')
        self.out_mut = tf.keras.layers.Dense(MUT_DIM, activation='linear')
        self.out_psi = tf.keras.layers.Dense(PSI_DIM, activation='linear')
        self.out_expr_sec = tf.keras.layers.Dense(EXPR_DIM, activation='linear')
        self.out_mut_sec = tf.keras.layers.Dense(MUT_DIM, activation='linear')
        self.out_psi_sec = tf.keras.layers.Dense(PSI_DIM, activation='linear')
        
                
    def call(self, inputs):
        
        out_expr = self.out_expr(inputs[0])
        out_mut = self.out_mut(inputs[0])
        out_psi = self.out_psi(inputs[0])
        out_expr_sec = self.out_expr_sec(inputs[1])
        out_mut_sec = self.out_mut_sec(inputs[2])
        out_psi_sec = self.out_psi_sec(inputs[3])
        d_expr = tf.math.add(out_expr, out_expr_sec)
        d_mut = tf.math.sigmoid(tf.math.add(out_mut, out_mut_sec))
        d_psi = tf.math.add(out_psi, out_psi_sec)
        
                
        return d_expr, d_mut, d_psi

class ModifyKLWeight(tf.keras.callbacks.Callback):
    
    def __init__(self):
        super(ModifyKLWeight, self).__init__()
        
    def on_epoch_end(self, epoch, logs=None):
        if not hasattr(self.model, 'custom_weight'):
            raise ValueError('Does not have weight attribute for KL')

        custom_weight = tf.keras.backend.get_value(self.model.custom_weight)
        if (epoch % 5 == 0) & (epoch > 1):
            tf.keras.backend.set_value(self.model.custom_weight, 0.0)
        else:
            if(custom_weight < 1.0): # Keep in case threshold needed
                tf.keras.backend.set_value(self.model.custom_weight, tf.keras.backend.get_value(self.model.custom_weight) + 0.1)
        print("\nEpoch %d: KL Weight is %.4f." % (epoch, tf.keras.backend.get_value(self.model.custom_weight)))



def sp_ae(data=None, dim=None, doCV=False, n_epoch=None, m_save=False):
    
    # Autoencoder
    BATCH_SIZE = 64
    if n_epoch is None:
        NUM_EPOCHS = 10000
    else:
        NUM_EPOCHS = n_epoch
        
    LAT_DIM = dim
    N_DIM_EXP=data[0].shape[1]
    N_DIM_MUT=data[1].shape[1]
    N_DIM_PSI=data[2].shape[1]
            
    # Model
    enc = EncoderModel(LAT_DIM)
    dec = DecoderModel(N_DIM_EXP, N_DIM_MUT, N_DIM_PSI)
    vae = VariationalAutoEncoder(enc, dec)
    
    if doCV is True:
        # Callback
        red_lr = tf.keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.1,patience=35,verbose=1,mode='auto', min_lr=1e-12)
        cust_cback = ModifyKLWeight()

        # Train
        rand_idx = np.random.permutation(data[0].shape[0])
        optimizer = tf.keras.optimizers.Adam(learning_rate=1e-4)
        vae.compile(optimizer=optimizer, loss=[cust_loss_mse,cust_loss,cust_loss_mse])
        hist = vae.fit([data[0][rand_idx,], data[1][rand_idx,], data[2][rand_idx,]], [data[0][rand_idx,], data[1][rand_idx,], data[2][rand_idx,]], epochs=NUM_EPOCHS, batch_size=BATCH_SIZE, shuffle=True, validation_split=0.2, callbacks=[cust_cback, red_lr])
    else:
        cust_cback = ModifyKLWeight()

        rand_idx = np.random.permutation(data[0].shape[0])
        optimizer = tf.keras.optimizers.Adam(learning_rate=1e-4)
        vae.compile(optimizer=optimizer, loss=[cust_loss_mse,cust_loss,cust_loss_mse])
        hist = vae.fit([data[0][rand_idx,], data[1][rand_idx,], data[2][rand_idx,]], [data[0][rand_idx,], data[1][rand_idx,], data[2][rand_idx,]], epochs=NUM_EPOCHS, batch_size=BATCH_SIZE, shuffle=True, callbacks=[cust_cback])
        
    # Save
    if m_save:
        enc.save('/data/splicing_enc_model')
        dec.save('/data/splicing_dec_model')

    # Clear
    tf.keras.backend.clear_session()
    
    return hist.history



if __name__ == '__main__':    
    
    l_wd = 'Research/splicing'
    
    # Load Data
    psi = scale(np.loadtxt('{}/data/psi_mat.tsv.gz'.format(l_wd), delimiter='\t', dtype=np.float64))
    gexp = scale(np.loadtxt('{}/data/expr_mat.tsv.gz'.format(l_wd), delimiter='\t', dtype=np.float64))
    mut = np.loadtxt('{}/data/mut_mat.tsv.gz'.format(l_wd), delimiter='\t', dtype=np.float64)

    # AE
    test_dim = np.asarray([16, 32, 64, 128, 256, 512])
    cv_res_psi = []
    cv_res_expr = []
    cv_res_mut = []
    for dim in test_dim:
        train_hist = sp_ae(data=[gexp, mut, psi], dim=dim, doCV=True, n_epoch=10)
        cv_res_expr.append(train_hist['val_output_1_loss'])
        cv_res_mut.append(train_hist['val_output_2_loss'])
        cv_res_psi.append(train_hist['val_output_3_loss'])
        
    #np.save('{}/data/CV_PSI_VAE.npy'.format(l_wd), arr=cv_res_psi)
    #np.save('{}/data/CV_EXPR_VAE.npy'.format(l_wd), arr=cv_res_expr)
    #np.save('{}/data/CV_MUT_VAE.npy'.format(l_wd), arr=cv_res_mut)
