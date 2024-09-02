// 总结内容： confirmSteps 绑定js函数 initOptionsUI
// // initOptionsUI 绑定js函数 init_type_div 绑定路由函数 /init_type_div 绑定函数 get_init_div(stepId, step_type, step_params, options_list, mode)
// // initOptionsUI 绑定js函数 update_my_options 绑定路由函数 /update_my_options 绑定函数 get_result_option(step_id, step_type, select_option,mode,temp_options)
// // initOptionsUI 绑定js函数 index_select 绑定路由函数 /index_select 绑定函数 get_type_index(stepId,step_type)


document.getElementById('confirmBtn').addEventListener('click', confirmSteps);
async function confirmSteps() {
    confirmBtn.disabled = true;  // 禁用确认按钮
    for (let i = 0; i < steps.length; i++) {steps[i].id = i;}
    for (const step of steps) { await initOptionsUI(step.id); } // 初始每一行的选项卡
    $('#selected-steps').sortable('disable'); // 禁用左侧按钮卡排序
    document.getElementById('submitBtn').disabled = false; // 允许提交
}


async function initOptionsUI(stepId) {
    const step = steps.find(step => step.id === stepId);
    const temp_data = await init_type_div(stepId,temp_options) // 获取给定step的输入框+每个选项卡的div内容，获取给定step的推荐选项字典
    document.getElementById('inputArea').insertAdjacentHTML('beforeend', temp_data.div_string); //添加给定step的div
    step.params = document.querySelector(`#input-step-${stepId} input`).value;      // 初始化step.params

    for (let key in temp_data.rec_option) {
        const selectElement = document.querySelector(`#input-step-${stepId} select#select-${stepId}-${key}`)
        const index = Array.from(selectElement.options).findIndex(option => option.value === temp_data.rec_option[key]);
        if (index !== -1) { selectElement.selectedIndex = index; } // 初始化展示的选项
        step.select_option[key] = selectElement.value; }                                    // 初始化step.select_option
    temp_options = await update_my_options(stepId);  //更新提供下一行的总选项

    if (step.type.slice(0, 3) === 'map') {
        const index_options =  await index_select(stepId);
        document.getElementById('inputArea').insertAdjacentHTML('beforeend', index_options);
        const selectElement = document.querySelector(`#index-step-${stepId} select`)
        selectElement.onchange = () => { step.select_index=selectElement.value; };
        step.select_index = selectElement.value;                                           // 初始化step.select_index
    }
}


async function init_type_div(stepId,options) {
    const step = steps.find(s => s.id === stepId);
    try {
        const response = await fetch('/init_type_div', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ step_id:stepId ,step_type: step.type, step_params: step.params,mode:step.mode,options:options })
        });
        return (await response.json());
    } catch (error) { console.error('Error:', error); }
}
async function update_my_options(stepId) {
    const step = steps.find(s => s.id === stepId);
    try {
        const response = await fetch('/update_my_options', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ step_id: step.id, step_type: step.type, select_option: step.select_option, mode: step.mode, temp_options: temp_options })
        });
        return  (await response.json()).temp_options;
    } catch (error) { console.error('Error:', error); }
}
async function index_select(stepId) {
    const step = steps.find(s => s.id === stepId);
    try {
        const response = await fetch('/index_select', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ stepId: stepId ,step_type: step.type })
        });
        return (await response.json()).index_options;
    } catch (error) { console.error('Error:', error); }
}

