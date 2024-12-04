

const presets_dict = {
    "preset_WGS_pair": "H4sIAOH7MWcC_7WQQQ7CIBBFr0JYm6agK3e6d-eONGRaoRKLhWlrEw13dzSRXqAuhsn_eeQFlCo3jJ_ReTbSoVvoejScus8EcEj7xW0UfM8UR5i1KGws2ievCLBR5l7--rTcKDNTJuKVoOwhsHqG1RwiM-LrkItDD8HdjLv_zbWlfDycGBrfP4y-TKFzDYxmWMFYg8__tiOyH68GWZygc_Q83ZAQV9PIVFVvMbmh9wwCAAA=",
    "preset_WGS_single": "H4sIADv8MWcC_52PsQ6CMBCGX-XSmRBAJzfd3dyapjmw1MaWllIk0fDuHg64isNdcn_yfbmf8yIDdonGQaIlNVofFaNsmcF02i7Xi7U9OwBnEae87XP9ZGJe42IWGfCSAIcB6gk38uWHr768HIK5K9P95dkRcDqeISrnH0pex2BNg0kNv9tqdGutPRE-3VSEfkRr6EPZkDdutVWzEG88ZkVXbgEAAA=="
};











async function initialPreset(base64_encoded_data) {
    try {
        const responseData = await postData('/batch_code/get_preset', {base64_encoded_data});
        responseData.forEach(element => {addStepOnly(element.type)});
        steps = responseData
        initOptionsUI(responseData)
    } catch (error) { console.error('Error initialPreset:', error);}
    disableAllButtons()
    confirmBtn.disabled = true;
    $('#selected-steps').sortable('disable');
    document.getElementById('submitBtn').disabled = false;
}

